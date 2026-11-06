# =============================================================================
# 04_table1.R
# Table 1. Characteristics of the 249 patients, overall and by current endemic status
#
# Run from the project root (the folder that contains data/ and scripts/).
# Input : data/clean/ipd_clean.csv
# Output: outputs/tables/Table1.docx   (Word table for the manuscript)
#         outputs/tables/Table1.csv    (same content as plain text)
# Packages: gtsummary, flextable, smd
#
# Descriptive only: median (IQR) and n (%). Groups are compared with standardised
# mean differences (SMD), not P values. Missing values are shown as their own row.
# =============================================================================

library(gtsummary)
library(flextable)

ipd <- read.csv(file.path("data", "clean", "ipd_clean.csv"),
                na.strings = "NA", stringsAsFactors = FALSE)
stopifnot(nrow(ipd) == 249)

# -----------------------------------------------------------------------------
# 1. Labelled factors (coding as in the data dictionary)
# -----------------------------------------------------------------------------
fac <- function(x, codes, labels, na_label = NULL) {
  f <- factor(x, levels = codes, labels = labels)
  if (!is.null(na_label)) {               # make missing an explicit category
    f <- addNA(f); levels(f)[is.na(levels(f))] <- na_label
  }
  f
}

# Diagnoses before confirmation: one group per patient, by priority
#   alternative diagnosis > brucellosis treated ineffectively > no diagnosis named >
#   brucellosis as the first diagnosis; "U" (not reported) becomes missing
PREDX_LEVELS <- c("Alternative diagnosis",
                  "Brucellosis diagnosed, treatment ineffective",
                  "No diagnosis named",
                  "Brucellosis was the first diagnosis")
predx_group <- function(x) {
  codes <- strsplit(ifelse(is.na(x), "", x), ";")
  g <- vapply(codes, function(cc) {
    cc <- trimws(cc)
    if (length(cc) == 0 || identical(cc, "") || identical(cc, "U")) return(NA_character_)
    if (any(cc %in% c("D", "T", "P", "N", "R", "E", "L"))) return(PREDX_LEVELS[1])
    if ("B1" %in% cc) return(PREDX_LEVELS[2])
    if (any(cc %in% c("S", "V"))) return(PREDX_LEVELS[3])
    if (identical(cc, "B0")) return(PREDX_LEVELS[4])
    stop("Unexpected PreDx_Category: ", paste(cc, collapse = ";"))
  }, character(1))
  factor(g, levels = PREDX_LEVELS)
}

t1 <- data.frame(
  endemic   = factor(ipd$Endemic_current, levels = c(1, 0),
                     labels = c("Endemic", "Non-endemic")),
  age       = ipd$Age,
  sex       = fac(ipd$Sex, c(1, 0), c("Male", "Female")),
  pub_era   = fac(ipd$Imaging_Era, c(0, 1), c("Before 2000", "2000 or later")),
  report    = fac(ipd$Report_Type, c(1, 2), c("Case report", "Case series")),
  exposure  = fac(ipd$Exposure_Hx, c(1, 0), c("Yes", "No or not mentioned")),
  fever     = fac(ipd$Fever_Status, c(1, 0), c("Yes", "No or not mentioned"),
                  na_label = "Reported only for the whole series"),
  temp      = ipd$Fever_Exact,
  delay     = ipd$Delay_to_Dx,
  dept      = fac(ipd$Dept_Bias, 1:3,
                  c("Internal medicine or infectious diseases",
                    "Orthopaedics or neurosurgery",
                    "Other departments or multidisciplinary team")),
  region    = fac(ipd$Spine_Region, c(3, 4), c("Lumbar or lumbosacral only",
                                              "Lumbar and other spinal regions")),
  lesions   = cut(ipd$Lesion_Count, breaks = c(-Inf, 0, 1, 2, Inf),
                  labels = c("0", "1", "2", "3 or more")),
  epidural  = fac(ipd$Epidural_Abscess, c(1, 0), c("Yes", "No")),
  psoas     = fac(ipd$Psoas_Abscess, c(1, 0), c("Yes", "No")),
  neuro     = fac(ipd$Neuro_Deficit, c(1, 0), c("Yes", "No")),
  extraspin = fac(ipd$Extra_Spinal_Comp, c(1, 0), c("Yes", "No")),
  crp       = ipd$CRP_Base,
  esr       = ipd$ESR_Base,
  regimen   = fac(ipd$Medical_Regimen, 1:3,
                  c("Dual therapy", "Aminoglycoside-containing", "Quinolone-intensified")),
  surgery   = fac(ipd$End_Status, c(1, 0, 2),
                  c("Spinal surgery", "Conservative treatment only",
                    "Died before surgery")),
  jbi       = ipd$JBI_Score,
  # Descriptive variables added 2026-10-21 (full breakdown in the supplementary table)
  prior     = fac(ipd$Prior_Brucellosis, c(1, 0), c("Yes", "No or not mentioned")),
  predx     = predx_group(ipd$PreDx_Category)
)

# Any code outside the expected values would silently become NA above; stop instead
check <- list(Spine_Region = c(3, 4), Medical_Regimen = 1:3, End_Status = 0:2, Dept_Bias = 1:3)
for (v in names(check)) {
  bad <- !is.na(ipd[[v]]) & !(ipd[[v]] %in% check[[v]])
  if (any(bad)) stop(v, " has unexpected codes: ", paste(unique(ipd[[v]][bad]), collapse = ", "))
}
lc <- ipd$Lesion_Count
bad <- !is.na(lc) & (lc < 0 | lc != round(lc))
if (any(bad)) stop("Lesion_Count must be a whole number >= 0: ", paste(unique(lc[bad]), collapse = ", "))
cat("Involved vertebrae = 0 (disc or epidural space only):",
    paste(ipd$Patient_ID[which(lc == 0)], collapse = ", "), "\n")

labels <- list(
  age = "Age, years", sex = "Sex", pub_era = "Year of publication",
  report = "Type of report", exposure = "Animal or raw dairy exposure",
  fever = "Fever during the current illness", temp = "Temperature at admission, °C",
  delay = "Diagnostic delay, months", dept = "Reporting department (corresponding author)",
  region = "Spinal involvement", lesions = "Number of involved vertebrae",
  epidural = "Epidural abscess", psoas = "Psoas abscess",
  neuro = "Neurological deficit", extraspin = "Severe extra-spinal complication",
  crp = "C-reactive protein, mg/L", esr = "Erythrocyte sedimentation rate, mm/h",
  regimen = "Antibiotic regimen", surgery = "Treatment", jbi = "JBI score (0–8)",
  prior = "Brucellosis diagnosed before the onset of back pain",
  predx = "Diagnosis before confirmation"
)

# -----------------------------------------------------------------------------
# 2. Table
# -----------------------------------------------------------------------------
# Absolute SMD, two decimals
smd_fmt <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.2f", abs(x)))

tbl <- tbl_summary(
  t1, by = endemic, label = labels,
  type = list(jbi ~ "continuous"),
  statistic = list(all_continuous() ~ "{median} ({p25}–{p75})",
                   all_categorical() ~ "{n} ({p}%)"),
  # whole numbers for age and JBI score (as recorded); 1 decimal for the rest
  digits = list(all_continuous() ~ 1, c(age, jbi) ~ 0, all_categorical() ~ c(0, 1)),
  missing = "ifany", missing_text = "Not reported"
) |>
  add_overall(last = FALSE) |>
  add_difference(test = everything() ~ "smd", estimate_fun = everything() ~ smd_fmt) |>
  modify_column_hide(columns = any_of(c("conf.low", "conf.high", "p.value"))) |>
  modify_header(label = "**Characteristic**", estimate = "**SMD**") |>
  modify_footnote(all_stat_cols() ~ paste(
    "Median (interquartile range) or n (%).",
    "No or not mentioned: history and symptom variables (exposure, fever and previous",
    "brucellosis) were coded as absent when the report did not mention them.",
    "Not reported: the report gave no information for the patient, or gave it only for the",
    "whole case series; shown as a separate row, excluded from the denominators of percentages",
    "and, for variables used in a model, from that model (complete-case analysis). Fever",
    "reported only for the whole series is shown as a separate category with all patients as the",
    "denominator, and is excluded from the main model.",
    "Involved vertebrae = 0: infection confined to the disc space or epidural space,",
    "without vertebral body involvement.",
    "Temperature at admission: highest value if several were reported.",
    "Diagnosis before confirmation: diagnoses recorded between the onset of back pain and",
    "confirmation of brucellosis with start of effective treatment; one group per patient, in",
    "the order alternative diagnosis, brucellosis treated ineffectively, no diagnosis named",
    "(symptomatic or empirical treatment only, or seen by a doctor with no diagnosis reported),",
    "brucellosis as the first diagnosis; percentages are among patients whose report described",
    "this period. Table S4 lists all categories, counting each patient in every category that applies.")) |>
  modify_footnote(estimate ~ paste(
    "Absolute standardised mean difference between endemic and non-endemic groups;",
    "values above 0.1 suggest imbalance.")) |>
  bold_labels()

# The SMD method adds a 'CI = Confidence Interval' abbreviation; no CI is shown, so drop it
tbl <- tryCatch(remove_abbreviation(tbl, "CI = Confidence Interval"),
                error = function(e) tbl)

print(tbl)

# -----------------------------------------------------------------------------
# 3. Save
# -----------------------------------------------------------------------------
out_dir <- file.path("outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(as_tibble(tbl, col_labels = TRUE)),
          file.path(out_dir, "Table1.csv"), row.names = FALSE)
ft <- as_flex_table(tbl) |>
  set_caption(caption = paste(
    "Table 1. Characteristics of", nrow(ipd), "patients with lumbar brucellar spondylitis,",
    "overall and by current endemic status")) |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 9, part = "all") |>
  padding(padding.top = 1, padding.bottom = 1, part = "body") |>   # tighter rows
  line_spacing(space = 1, part = "all") |>
  autofit()
save_as_docx(ft, path = file.path(out_dir, "Table1.docx"))
cat("Saved Table 1 to", normalizePath(out_dir), "\n")
