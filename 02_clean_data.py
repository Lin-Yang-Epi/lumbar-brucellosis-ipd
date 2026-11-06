"""
02_clean_data.py
Clean the lumbar brucellar spondylitis IPD file and write analysis-ready CSVs for R.

Input : data/raw/data_cleaned_final.xlsx   (sheet "Sheet1", 249 patients; sheet "Excluded_cases" lists
                                             patients removed at case level, not read here)
Output: data/clean/ipd_clean.csv           analysis file (NA written as "NA")
        data/clean/ipd_dictionary.csv      variable dictionary for the clean file
        data/clean/ipd_missingness.csv     missing count per variable
        data/clean/ipd_checks.txt          validation report

Rules
  * Empty or whitespace-only cells, "NA", "N/A" and "." become true missing values.
  * No imputation. Every value is exactly as reported in the source publication.
  * Raw columns are kept unchanged. Derived columns are added with new names.
  * SAT_Titer, Coombs_Test and Gold_Standard stay in the raw file only.

Run from the project root:  python scripts/python/02_clean_data.py
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
RAW_FILE = ROOT / "data" / "raw" / "data_cleaned_final.xlsx"
OUT_DIR = ROOT / "data" / "clean"
ENDEMIC_FILE = ROOT / "data" / "raw" / "endemic_classification.csv"

EXPECTED_N = 249
IMAGING_CUTOFF_YEAR = 2000      # MRI / CT widely available from this year
ZERO_TIME_REPLACEMENT = 0.1     # months; Time_to_Event = 0 means immediate surgery

NA_STRINGS = {"", "NA", "N/A", "NAN", "NULL", "NONE", "."}

# Columns not carried into the analysis file (diagnostic information only)
DROP_FROM_ANALYSIS = ["SAT_Titer", "Coombs_Test", "Gold_Standard"]

# Columns kept in the raw file only, deliberately not used (no warning)
IGNORED_COLUMNS = ["首诊诊断", "备注"]   # free-text notes, kept in the raw file only

# Header variants that should map to the dictionary name
RENAME = {"Coombs test": "Coombs_Test", "coombs test": "Coombs_Test"}

# --------------------------------------------------------------------------
# Data dictionary: name -> (type, allowed values or range, label)
#   type: text | int | cat | num
# --------------------------------------------------------------------------
SPEC = {
    "Study_ID":          ("text", None, "Study unique ID"),
    "Patient_ID":        ("text", None, "Patient unique ID"),
    "Pub_Year":          ("int",  (1900, 2026), "Year of publication"),
    "Region":            ("text", None, "Country of diagnosis"),
    "Dept_Bias":         ("cat",  {1: "Internal medicine / infectious diseases",
                                   2: "Orthopaedics / neurosurgery",
                                   3: "Other departments / multidisciplinary team"},
                          "Department of the corresponding author (reporting department; descriptive only)"),
    "Age":               ("num",  (0, 110), "Age, years"),
    "Sex":               ("cat",  {0: "Female", 1: "Male"}, "Sex"),
    "Exposure_Hx":       ("cat",  {0: "Not mentioned or none", 1: "Clear exposure"},
                          "Brucella exposure history"),
    "Delay_to_Dx":       ("num",  (0, 240), "Onset of back pain (first lumbar-related symptom of the "
                                             "current illness) to diagnosis and start of effective "
                                             "treatment, months"),
    "Fever_Exact":       ("num",  (35, 42), "Admission temperature, deg C (highest if several)"),
    "Fever_Status":      ("cat",  {0: "No or not mentioned", 1: "Yes"},
                          "Fever during the current illness (NA = reported only as a group proportion)"),
    "Neuro_Deficit":     ("cat",  {0: "No", 1: "Yes"},
                          "Neurological deficit, strict definition (muscle weakness or cauda equina/"
                          "conus syndrome on examination at diagnosis; main definition)"),
    "Neuro_Deficit_Broad": ("cat", {0: "No", 1: "Yes"},
                          "Neurological deficit, broad definition (strict definition or any objective "
                          "neurological sign: segmental atrophy, reduced/absent reflexes, upper motor "
                          "neuron signs, objective sensory loss; unplanned analysis)"),
    "Extra_Spinal_Comp": ("cat",  {0: "Isolated spine infection",
                                   1: "Severe systemic complication"},
                          "Extra-spinal complication"),
    "Spine_Region":      ("cat",  {1: "Cervical", 2: "Thoracic", 3: "Lumbar / sacral",
                                   4: "Multiregional / skip"}, "Spinal region"),
    "Lesion_Count":      ("int",  (0, 30), "Number of involved vertebral segments"),
    "Epidural_Abscess":  ("cat",  {0: "No", 1: "Yes"}, "Epidural abscess"),
    "Psoas_Abscess":     ("cat",  {0: "No", 1: "Yes"}, "Psoas abscess"),
    "CRP_Base":          ("num",  (0, 2000), "Baseline CRP, mg/L (descriptive only)"),
    "ESR_Base":          ("num",  (0, 200), "Baseline ESR, mm/h (descriptive only)"),
    "SAT_Titer":         ("num",  (0, None), "Serum agglutination titre (raw file only)"),
    "Coombs_Test":       ("num",  (0, None), "Coombs titre (raw file only)"),
    "Gold_Standard":     ("cat",  {1: "Serology only", 2: "Blood culture positive",
                                   3: "Focal biopsy culture / mNGS positive"},
                          "Diagnostic basis (raw file only)"),
    "Medical_Regimen":   ("cat",  {1: "Dual therapy", 2: "Aminoglycoside-containing",
                                   3: "Quinolone-intensified"}, "Medical regimen"),
    "Time_to_Event":     ("num",  (0, 240), "Months from start of conservative treatment "
                                             "to surgery / death / last follow-up"),
    "End_Status":        ("cat",  {0: "Censored", 1: "Surgery",
                                   2: "Competing event (pre-operative death)"},
                          "Outcome status"),
    "Surg_Type":         ("cat",  {1: "Abscess drainage / decompression",
                                   2: "Internal fixation / fusion"},
                          "Surgery type (NA = no surgery)"),
    # Optional columns: the script runs without them, and uses them once they are added
    "Report_ID":         ("text", None, "Publication label. Leave blank to derive it from Patient_ID; "
                                        "fill only where the Patient_ID letter does not identify the paper"),
    "Report_N":          ("int",  (1, 1000), "Number of spinal brucellosis patients reported in the "
                                             "original publication (eligibility: no more than 50)"),
    "Report_N_total":    ("int",  (1, 100000), "Total number of patients reported in the original "
                                               "publication, any site (defines Report_Type)"),
    "JBI_Score":        ("int",  (0, 8), "Number of the 8 JBI items answered Yes"),
    "JBI_Overall":       ("text", None, "JBI overall decision (Include / Exclude / Seek further info)"),
    # Descriptive only, not in any model (added 2026-10-21)
    "Prior_Brucellosis": ("cat",  {0: "No or not mentioned",
                                   1: "Brucellosis diagnosed before the onset of back pain"},
                          "Previous brucellosis (descriptive only)"),
    "PreDx_Reported":    ("cat",  {0: "Not reported", 1: "Reported"},
                          "Diagnoses and treatment before confirmation described in the report "
                          "(descriptive only)"),
    "PreDx_Category":    ("text", None, "Diagnoses before confirmation, codes separated by ';' "
                                        "(see PREDX_CODES; descriptive only)"),
}
OPTIONAL_COLS = ["Report_ID", "Report_N", "Report_N_total", "JBI_Score", "JBI_Overall",
                 "Neuro_Deficit_Broad"]
MAX_REPORT_N = 50   # eligibility: no more than 50 spinal brucellosis patients per publication
JBI_OVERALL_VALUES = {"include", "exclude", "seek further info"}

# Diagnoses before confirmation (PreDx_Category); one patient can have several
PREDX_CODES = {
    "D":  "Degenerative or other non-infectious spinal disease",
    "T":  "Spinal tuberculosis (suspected or treated)",
    "P":  "Pyogenic or unspecified spinal infection",
    "N":  "Tumour, metastasis or space-occupying lesion",
    "R":  "Rheumatic or immune-mediated disease",
    "E":  "Extra-spinal disease",
    "S":  "Symptomatic or empirical treatment only, no diagnosis named",
    "V":  "Seen by a doctor, diagnosis and treatment not reported",
    "L":  "Laboratory misidentification",
    "B1": "Brucellosis diagnosed but treatment ineffective",
    "B0": "First diagnosis was brucellosis",
    "U":  "Not reported",
}

DERIVED = {
    "Time_to_Event_adj": ("num", "Time_to_Event with 0 replaced by 0.1 months"),
    "Log_Delay":         ("num", "Natural log of Delay_to_Dx"),
    "Imaging_Era":       ("cat", f"0 = published before {IMAGING_CUTOFF_YEAR}; "
                                 f"1 = published {IMAGING_CUTOFF_YEAR} or later"),
    "Endemic_current":   ("cat", "Country status today: 1 = endemic; 0 = non-endemic "
                                 "(official brucellosis-free status or only sporadic cases)"),
    "Endemic_Pappas":    ("cat", "Country status per Pappas 2006: 1 = endemic; 0 = non-endemic; "
                                 "countries without Pappas data take the current status"),
    "Cluster_ID":        ("text", "Source publication. Report_ID where filled; otherwise Patient_ID "
                                  "with the trailing letter dropped (Farid_1963a -> Farid_1963)"),
    "Cluster_N":         ("int", "Number of patients from the same report within the analysis cohort"),
    "Report_Type":       ("cat", "Primary definition from Report_N_total: 1 = case report "
                                 "(1 patient in the publication); 2 = case series (2 or more)"),
    "Report_Type_alt":   ("cat", "Robustness definition from Report_N_total: 1 = 1 or 2 patients; "
                                 "2 = 3 or more patients"),
}

errors, warnings, notes = [], [], []


def err(msg):
    errors.append(msg)


def warn(msg):
    warnings.append(msg)


# --------------------------------------------------------------------------
# 1. Read
# --------------------------------------------------------------------------
if not RAW_FILE.exists():
    sys.exit(f"Input file not found: {RAW_FILE}")

df = pd.read_excel(RAW_FILE, sheet_name="Sheet1", dtype=object)
df.columns = [RENAME.get(str(c).strip(), str(c).strip()) for c in df.columns]

missing_cols = [c for c in SPEC if c not in df.columns and c not in OPTIONAL_COLS]
missing_optional = [c for c in OPTIONAL_COLS if c not in df.columns]
extra_cols = [c for c in df.columns if c not in SPEC]
if missing_cols:
    err(f"Columns in dictionary but not in file: {missing_cols}")
if missing_optional:
    notes.append(f"Optional columns not in the file yet: {missing_optional}")
extra_cols = [c for c in extra_cols if c not in IGNORED_COLUMNS]
if extra_cols:
    warn(f"Columns in file but not in dictionary (ignored): {extra_cols}")
df = df[[c for c in SPEC if c in df.columns]].copy()

# --------------------------------------------------------------------------
# 2. Blank cells -> true NA
# --------------------------------------------------------------------------
def to_na(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return pd.NA
    if isinstance(x, str):
        s = x.strip()
        return pd.NA if s.upper() in NA_STRINGS else s
    return x


df = df.apply(lambda col: col.map(to_na))

# --------------------------------------------------------------------------
# 3. Types
# --------------------------------------------------------------------------
for col, (kind, rule, _label) in SPEC.items():
    if col not in df.columns:
        continue
    if kind == "text":
        s = df[col].astype("string").str.strip()
        # numeric IDs read from Excel as 1.0 -> "1"
        s = s.str.replace(r"\.0$", "", regex=True)
        df[col] = s
    else:
        num = pd.to_numeric(df[col], errors="coerce")
        bad = df[col].notna() & num.isna()
        if bad.any():
            err(f"{col}: non-numeric entries at rows "
                f"{(df.index[bad] + 2).tolist()} -> {df.loc[bad, col].tolist()}")
        if kind in ("int", "cat"):
            frac = num.notna() & (num != num.round())
            if frac.any():
                err(f"{col}: non-integer values at rows {(df.index[frac] + 2).tolist()} "
                    f"-> {num[frac].tolist()}")
            # set non-integers to NA so the run can finish and list every error
            df[col] = num.mask(frac).astype("Int64")
        else:
            df[col] = num.astype("Float64")

# --------------------------------------------------------------------------
# 4. Validation
# --------------------------------------------------------------------------
# 4a. size and IDs
if len(df) != EXPECTED_N:
    err(f"Expected {EXPECTED_N} rows, found {len(df)}")
if df["Patient_ID"].isna().any():
    err("Patient_ID is missing for some rows")
dup = df["Patient_ID"][df["Patient_ID"].duplicated(keep=False) & df["Patient_ID"].notna()]
if len(dup):
    warn(f"Patient_ID is not unique ({dup.nunique()} IDs repeated): "
         f"{sorted(dup.unique().tolist())[:15]}")
    # Show the repeated rows (Excel row numbers) so they can be traced to the source
    show = df.loc[dup.index, ["Patient_ID", "Study_ID", "Age", "Sex", "Delay_to_Dx"]].copy()
    show.insert(0, "excel_row", show.index + 2)
    show = show.sort_values(["Patient_ID", "excel_row"])
    notes.append("Rows sharing a Patient_ID:\n" + show.to_string(index=False))
if df["Study_ID"].isna().any():
    err("Study_ID is missing for some rows")
notes.append(f"Rows: {len(df)}; countries: {df['Region'].nunique()}")

# 4b. allowed values and ranges
for col, (kind, rule, _label) in SPEC.items():
    if col not in df.columns or rule is None:
        continue
    s = df[col]
    if kind == "cat":
        bad = s.notna() & ~s.isin(list(rule.keys()))
        if bad.any():
            err(f"{col}: values outside {sorted(rule)} at rows "
                f"{(df.index[bad] + 2).tolist()} -> {s[bad].tolist()}")
    else:
        lo, hi = rule
        if lo is not None and (s < lo).any():
            err(f"{col}: values below {lo} at rows {(df.index[s < lo] + 2).tolist()}")
        if hi is not None and (s > hi).any():
            err(f"{col}: values above {hi} at rows {(df.index[s > hi] + 2).tolist()}")

# 4c. required variables without missing values
for col in ["Pub_Year", "Region", "Age", "Sex", "Delay_to_Dx", "Exposure_Hx"]:
    n = int(df[col].isna().sum())
    if n:
        err(f"{col}: {n} missing values in a variable that should be complete "
            f"(rows {(df.index[df[col].isna()] + 2).tolist()[:20]})")
# Variables allowed to be missing (group-level reporting, or not needed for the delay analysis)
allowed_na = {c: int(df[c].isna().sum()) for c in
              ["Fever_Status", "Epidural_Abscess", "Neuro_Deficit", "Neuro_Deficit_Broad", "Dept_Bias",
               "Time_to_Event", "End_Status"] if c in df.columns}
notes.append(f"Missing values allowed (count): {allowed_na}")
fs = {("NA" if pd.isna(k) else int(k)): int(v)
      for k, v in df["Fever_Status"].value_counts(dropna=False).sort_index().items()}
notes.append(f"Fever_Status (1 = yes, 0 = no/not mentioned, NA = group-level only): {fs}")
ex = {int(k): int(v) for k, v in df["Exposure_Hx"].value_counts().sort_index().items()}
notes.append(f"Exposure_Hx (1 = clear exposure, 0 = not mentioned or none): {ex}")
dp = {("NA" if pd.isna(k) else int(k)): int(v)
      for k, v in df["Dept_Bias"].value_counts(dropna=False).sort_index().items()}
notes.append(f"Dept_Bias (1 = internal medicine/ID, 2 = orthopaedics/neurosurgery, "
             f"3 = other/MDT): {dp}")

# 4d. logical consistency
def flag(mask):
    """Nullable comparisons give NA; treat NA as 'not flagged'."""
    return mask.fillna(False).astype(bool)


# Broad neurological definition must contain the strict one, with the same missing values
if "Neuro_Deficit_Broad" in df.columns:
    nd, nb = df["Neuro_Deficit"], df["Neuro_Deficit_Broad"]
    bad = nd.isna() != nb.isna()
    if bad.any():
        err(f"Neuro_Deficit and Neuro_Deficit_Broad differ in missingness at rows {(df.index[bad] + 2).tolist()}")
    bad = flag(nd == 1) & flag(nb != 1)
    if bad.any():
        err(f"Neuro_Deficit = 1 but Neuro_Deficit_Broad != 1 at rows {(df.index[bad] + 2).tolist()}")

surg = flag(df["End_Status"] == 1)
bad = surg & df["Surg_Type"].isna()
if bad.any():
    warn(f"Surgery recorded but Surg_Type missing at rows {(df.index[bad] + 2).tolist()}")
bad = flag(df["End_Status"] != 1) & df["Surg_Type"].notna()
if bad.any():
    warn(f"Surg_Type present without End_Status = 1 at rows {(df.index[bad] + 2).tolist()}")
bad = flag(df["Fever_Exact"].notna() & (df["Fever_Status"] == 0) & (df["Fever_Exact"] >= 37.5))
if bad.any():
    warn(f"Fever_Exact >= 37.5 but Fever_Status = 0 at rows {(df.index[bad] + 2).tolist()}")
bad = flag(df["Fever_Exact"].notna() & (df["Fever_Status"] == 1) & (df["Fever_Exact"] < 37.5))
if bad.any():
    warn(f"Fever_Status = 1 but Fever_Exact < 37.5 at rows {(df.index[bad] + 2).tolist()}")
bad = df["Fever_Exact"].notna() & df["Fever_Status"].isna()
if bad.any():
    warn(f"Fever_Exact present but Fever_Status missing at rows {(df.index[bad] + 2).tolist()}")
if flag(df["Time_to_Event"] == 0).any():
    notes.append(f"Time_to_Event = 0 in {int((df['Time_to_Event'] == 0).sum())} patients "
                 f"(immediate surgery); replaced by {ZERO_TIME_REPLACEMENT} in Time_to_Event_adj")
n_ce = int((df["End_Status"] == 2).sum())
counts = {int(k): int(v) for k, v in df["End_Status"].value_counts().sort_index().items()}
notes.append(f"End_Status counts: {counts}")
if n_ce != 2:
    warn(f"End_Status = 2 in {n_ce} patients; the re-adjudicated cohort should have 2 "
         f"(Hegde 2025, Sanchez-Gonzalez 1996)")
if df["Spine_Region"].isin([4]).any():
    notes.append("Spine_Region = 4 is 'multiregional / skip' in the dictionary "
                 "(the old Table 1 code labelled it 'Lumbosacral')")

# --------------------------------------------------------------------------
# 5. Derived variables (raw columns stay unchanged)
# --------------------------------------------------------------------------
df["Time_to_Event_adj"] = df["Time_to_Event"].mask(flag(df["Time_to_Event"] == 0),
                                                    ZERO_TIME_REPLACEMENT)
if (df["Delay_to_Dx"] <= 0).any():
    err("Delay_to_Dx <= 0 found; Log_Delay cannot be computed")
df["Log_Delay"] = np.log(df["Delay_to_Dx"].astype("float64")).astype("Float64")
df["Imaging_Era"] = (df["Pub_Year"] >= IMAGING_CUTOFF_YEAR).astype("Int64")

# Report-level cluster: patients from the same publication share a Patient_ID stem.
# Study_ID is one number per patient, so it cannot be used for clustering.
bad_id = ~df["Patient_ID"].str.contains(r"_\d{4}[a-z]*$", regex=True, na=False)
if bad_id.any():
    warn(f"Patient_ID not in the form Name_YYYY or Name_YYYYa at rows "
         f"{(df.index[bad_id] + 2).tolist()}")
stem = df["Patient_ID"].str.replace(r"(?<=\d{4})[a-z]+$", "", regex=True)
if "Report_ID" in df.columns:
    # Report_ID overrides the derived stem where it is filled
    df["Cluster_ID"] = df["Report_ID"].str.strip().where(df["Report_ID"].notna(), stem)
    n_over = int(df["Report_ID"].notna().sum())
    notes.append(f"Report_ID filled for {n_over} patients; the rest use the Patient_ID stem")
    # One Report_ID spread over different stems is suspicious
    spread = df[df["Report_ID"].notna()].groupby("Report_ID")["Patient_ID"].apply(
        lambda s: s.str.replace(r"(?<=\d{4})[a-z]+$", "", regex=True).nunique())
    if (spread > 1).any():
        warn(f"Report_ID used across different Patient_ID stems: {spread[spread > 1].index.tolist()}")
    # A Report_ID that duplicates another report's derived stem would merge two papers
    clash = set(df.loc[df["Report_ID"].notna(), "Report_ID"]) & set(stem[df["Report_ID"].isna()])
    if clash:
        warn(f"Report_ID equals the stem of other patients without a Report_ID: {sorted(clash)}")
else:
    df["Cluster_ID"] = stem
df["Cluster_N"] = df.groupby("Cluster_ID")["Patient_ID"].transform("size").astype("Int64")
sizes = df.drop_duplicates("Cluster_ID")["Cluster_N"].value_counts().sort_index()
notes.append(f"Reports (clusters): {df['Cluster_ID'].nunique()}; patients per report: "
             f"{ {int(k): int(v) for k, v in sizes.items()} }")
multi = df[df["Cluster_N"] > 1].assign(excel_row=lambda d: d.index + 2)
if len(multi):
    summ = (multi.groupby("Cluster_ID")
            .agg(n=("Patient_ID", "size"),
                 patient_ids=("Patient_ID", lambda s: ", ".join(s) if len(s) <= 6
                              else f"{s.iloc[0]} ... {s.iloc[-1]}"),
                 excel_rows=("excel_row", lambda r: ", ".join(map(str, r)) if len(r) <= 6
                             else f"{r.min()}-{r.max()}"))
            .reset_index())
    notes.append("Reports with more than one patient (check these are the same publication):\n"
                 + summ.to_string(index=False))
if df["Study_ID"].nunique() == len(df):
    notes.append("Study_ID is unique for every patient, so it is not a cluster variable; "
                 "use Cluster_ID")

# Report size (Report_N, spinal cases) and report type (Report_N_total, all cases)
def rows_of(mask):
    mask = mask.fillna(False).astype(bool)
    return (df.index[mask] + 2).tolist()


def check_constant_within_report(col):
    varying = df.groupby("Cluster_ID")[col].nunique()
    varying = varying[varying > 1].index.tolist()
    if varying:
        err(f"{col} differs between patients of the same Cluster_ID: {varying}. "
            f"If these are different papers, fill Report_ID for them (one label per paper)")


if "Report_N" in df.columns:
    if df["Report_N"].isna().any():
        err(f"Report_N missing at rows {rows_of(df['Report_N'].isna())}")
    # Eligibility: no more than 50 spinal brucellosis patients in the publication
    too_big = df["Report_N"] > MAX_REPORT_N
    if too_big.any():
        err(f"Report_N (spinal cases) > {MAX_REPORT_N}, fails the eligibility criterion, rows "
            f"{rows_of(too_big)}, Patient_ID {df.loc[too_big, 'Patient_ID'].tolist()}")
    check_constant_within_report("Report_N")
    small = df["Report_N"] < df["Cluster_N"]
    if small.any():
        err(f"Report_N is smaller than the number of patients taken from that Cluster_ID, "
            f"rows {rows_of(small)}. If these are different papers, "
            f"fill Report_ID for them (one label per paper)")

if "Report_N_total" in df.columns:
    if df["Report_N_total"].isna().any():
        err(f"Report_N_total missing at rows {rows_of(df['Report_N_total'].isna())}")
    check_constant_within_report("Report_N_total")
    if "Report_N" in df.columns:
        less = df["Report_N_total"] < df["Report_N"]
        if less.any():
            err(f"Report_N_total (all cases) is smaller than Report_N (spinal cases) "
                f"at rows {rows_of(less)}")
        differ = df.loc[df["Report_N_total"] > df["Report_N"], ["Cluster_ID", "Report_N",
                                                                "Report_N_total"]]
        if len(differ):
            differ = differ.drop_duplicates("Cluster_ID")
            notes.append("Reports that also include non-spinal patients "
                         "(Report_N_total > Report_N):\n" + differ.to_string(index=False))
    df["Report_Type"] = df["Report_N_total"].map(
        lambda x: pd.NA if pd.isna(x) else (1 if x == 1 else 2)).astype("Int64")
    df["Report_Type_alt"] = df["Report_N_total"].map(
        lambda x: pd.NA if pd.isna(x) else (1 if x <= 2 else 2)).astype("Int64")
    notes.append(f"Report_Type (1 vs 2+ patients): "
                 f"{int((df['Report_Type'] == 1).sum())} patients from case reports, "
                 f"{int((df['Report_Type'] == 2).sum())} from case series")
    notes.append(f"Report_Type_alt (1-2 vs 3+ patients): "
                 f"{int((df['Report_Type_alt'] == 1).sum())} patients from case reports, "
                 f"{int((df['Report_Type_alt'] == 2).sum())} from case series")
elif "Report_N" in df.columns:
    notes.append("Report_N_total not in the file yet, so Report_Type is not derived")

# JBI columns
if "JBI_Overall" in df.columns:
    bad = df["JBI_Overall"].notna() & ~df["JBI_Overall"].str.lower().isin(JBI_OVERALL_VALUES)
    if bad.any():
        err(f"JBI_Overall has unexpected values at rows {(df.index[bad] + 2).tolist()}")
    notes.append(f"JBI_Overall: {df['JBI_Overall'].value_counts(dropna=False).to_dict()}")
if "JBI_Score" in df.columns:
    js = {("NA" if pd.isna(k) else int(k)): int(v)
          for k, v in df["JBI_Score"].value_counts(dropna=False).sort_index().items()}
    notes.append(f"JBI_Score: {js}")

# Descriptive variables: previous brucellosis and diagnoses before confirmation
if "Prior_Brucellosis" in df.columns:
    notes.append(f"Prior_Brucellosis = 1: {int((df['Prior_Brucellosis'] == 1).sum())} patients")
if "PreDx_Category" in df.columns:
    codes = df["PreDx_Category"].fillna("").str.split(";")
    bad = codes.map(lambda xs: any(x.strip() not in PREDX_CODES for x in xs))
    if bad.any():
        err(f"PreDx_Category has unknown codes at rows {rows_of(bad)}")
    if "PreDx_Reported" in df.columns:
        is_u = df["PreDx_Category"].str.strip() == "U"
        mism = flag((df["PreDx_Reported"] == 0) != is_u)
        if mism.any():
            err(f"PreDx_Reported does not match PreDx_Category (U = not reported) at rows {rows_of(mism)}")
        n_rep = int((df["PreDx_Reported"] == 1).sum())
        cnt = {k: int(codes.map(lambda xs: k in [x.strip() for x in xs]).sum()) for k in PREDX_CODES}
        notes.append(f"PreDx_Reported = 1: {n_rep} patients; PreDx_Category counts (a patient can have "
                     f"several): {cnt}")

# Endemic classification by country (see endemic_classification.csv for sources)
if not ENDEMIC_FILE.exists():
    err(f"Endemic classification file not found: {ENDEMIC_FILE}")
else:
    ec = pd.read_csv(ENDEMIC_FILE, dtype=str)
    ec["Region"] = ec["Region"].str.strip()
    code = {"Endemic": 1, "Non-endemic": 0}
    cur = ec.set_index("Region")["Current_class_proposed"].map(code)
    pap = ec.set_index("Region")["Pappas2006_class"].map(code).fillna(cur)  # No data -> current
    unmatched = sorted(set(df["Region"].dropna()) - set(cur.index))
    if unmatched:
        err(f"Region values not in the endemic classification file: {unmatched}")
    df["Endemic_current"] = df["Region"].map(cur).astype("Int64")
    df["Endemic_Pappas"] = df["Region"].map(pap).astype("Int64")
    changed = df.loc[flag(df["Endemic_current"] != df["Endemic_Pappas"]), "Region"]
    notes.append(f"Endemic_current: {int((df['Endemic_current'] == 1).sum())} endemic, "
                 f"{int((df['Endemic_current'] == 0).sum())} non-endemic")
    notes.append(f"Patients classified differently by Pappas 2006: {len(changed)} "
                 f"{changed.value_counts().to_dict()}")
    # Patient counts recorded in the classification file should match the data
    if "N_patients" in ec.columns:
        expected = pd.to_numeric(ec.set_index("Region")["N_patients"], errors="coerce")
        actual = df["Region"].value_counts()
        diff = {r: (int(actual.get(r, 0)), None if pd.isna(expected.get(r)) else int(expected.get(r)))
                for r in set(actual.index) | set(expected.dropna().index)
                if int(actual.get(r, 0)) != (0 if pd.isna(expected.get(r)) else int(expected.get(r)))}
        if diff:
            warn(f"N_patients in the classification file differs from the data "
                 f"(country: (data, file)): {diff}")
    notes.append(f"Countries: {df['Region'].nunique()}; patients per country: "
                 f"{df['Region'].value_counts().to_dict()}")

# Number of patients available for each pre-specified model (complete cases)
MAIN_COVARS = ["Log_Delay", "Age", "Sex", "Exposure_Hx", "Endemic_current",
               "Imaging_Era", "Fever_Status"]
SECONDARY_COVARS = ["Log_Delay", "Age", "Sex", "Endemic_current", "Imaging_Era"]
if "Endemic_current" in df.columns:
    n_main = int(df[MAIN_COVARS].notna().all(axis=1).sum())
    n_reduced = int(df[[c for c in MAIN_COVARS if c != "Fever_Status"]].notna().all(axis=1).sum())
    notes.append(f"Complete cases: main model {n_main}; reduced model without fever (S6) {n_reduced}")
    for outcome in ["Epidural_Abscess", "Neuro_Deficit", "Neuro_Deficit_Broad"]:
        if outcome not in df.columns:
            continue
        cc = df[SECONDARY_COVARS + [outcome]].dropna()
        notes.append(f"Complete cases: secondary model for {outcome} {len(cc)} "
                     f"(events {int((cc[outcome] == 1).sum())})")
    if "JBI_Score" in df.columns:
        n_s4 = int((df[MAIN_COVARS].notna().all(axis=1) & flag(df["JBI_Score"] == 8)).sum())
        notes.append(f"Complete cases with JBI_Score = 8 (S4): {n_s4}")

# --------------------------------------------------------------------------
# 6. Analysis file and outputs
# --------------------------------------------------------------------------
analysis = df.drop(columns=[c for c in DROP_FROM_ANALYSIS if c in df.columns])

OUT_DIR.mkdir(parents=True, exist_ok=True)
analysis.to_csv(OUT_DIR / "ipd_clean.csv", index=False, na_rep="NA", encoding="utf-8")

rows = []
for col in analysis.columns:
    if col in SPEC:
        kind, rule, label = SPEC[col]
        coding = ("; ".join(f"{k} = {v}" for k, v in rule.items())
                  if kind == "cat" else
                  (f"range {rule[0]} to {rule[1]}" if rule else ""))
        origin = "raw"
        if col == "PreDx_Category":
            coding = "; ".join(f"{k} = {v}" for k, v in PREDX_CODES.items())
    else:
        kind, label = DERIVED[col]
        coding, origin = "", "derived"
        if col == "Imaging_Era":
            coding = "0 = before cutoff year; 1 = cutoff year or later"
    rows.append({"variable": col, "type": kind, "label": label,
                 "coding": coding, "origin": origin,
                 "n_missing": int(analysis[col].isna().sum())})
pd.DataFrame(rows).to_csv(OUT_DIR / "ipd_dictionary.csv", index=False, encoding="utf-8")

miss = (analysis.isna().sum().rename("n_missing").to_frame()
        .assign(pct_missing=lambda d: (100 * d["n_missing"] / len(analysis)).round(1)))
miss.index.name = "variable"
miss.to_csv(OUT_DIR / "ipd_missingness.csv", encoding="utf-8")

report = ["VALIDATION REPORT", "=" * 60, ""]
report += ["NOTES"] + [f"  - {m}" for m in notes] + [""]
report += [f"WARNINGS ({len(warnings)})"] + [f"  - {m}" for m in warnings] + [""]
report += [f"ERRORS ({len(errors)})"] + [f"  - {m}" for m in errors] + [""]
report += ["MISSING VALUES (analysis file)", miss.to_string()]
(OUT_DIR / "ipd_checks.txt").write_text("\n".join(report), encoding="utf-8")

print("\n".join(report))
print(f"\nWritten to {OUT_DIR}")
if errors:
    sys.exit(1)
