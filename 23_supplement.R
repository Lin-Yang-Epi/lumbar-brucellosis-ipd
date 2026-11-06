# =============================================================================
# 23_supplement.R
# One Electronic Supplementary Material file: Tables S1-S21 and Figs. S1-S4,
# one item per page, A4 landscape, numbered in the order they are first cited
# in the paper (tables first, then figures).
#
#   Tables S1, S3, S5-S8  built here from their CSV files (caption + notes added
#                  here; each is also saved as its own .docx)
#   Other tables   read from the .rds file that save_docx_table() writes next
#                  to each .docx, so captions and notes are exactly those of the
#                  individual tables (scripts 07-15, 17-19, 21, 22 and 23)
#   Figs. S1-S4    the PNG files, with the legend below each figure
# The table of contents is added by hand in Word.
#
# Run from the project root, after scripts 02-19 and 21-23.
# Output: outputs/Supplementary_Material.docx
#         outputs/tables/TableS1_countries.docx, TableS3_JBI_items.docx,
#         TableS5_collinearity.docx ... TableS8_secondary_covariates.docx
# Packages: officer, flextable (plus those in 00_helpers.R)
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages({ library(officer); library(flextable) })

ipd  <- load_ipd()
cc   <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS)]), ]
n_cc <- nrow(cc); k_cc <- length(unique(cc$Cluster_ID))

rd <- function(f) {
  path <- file.path(OUT_TAB, f)
  if (!file.exists(path)) stop("Not found: ", path, " - run the script that makes it first")
  d <- read.csv(path, check.names = FALSE, colClasses = "character",
                na.strings = c("", "NA"), fileEncoding = "UTF-8")
  d[is.na(d)] <- ""
  d
}

# Stack groups as panels: a bold header row per group, the group's rows below
stack_panels <- function(d, group_col, keep) {
  parts <- list(); hdr <- integer(0); r <- 0
  for (g in unique(d[[group_col]])) {
    sub  <- d[d[[group_col]] == g, keep, drop = FALSE]
    head <- sub[1, , drop = FALSE]; head[1, ] <- ""; head[1, 1] <- g
    parts[[length(parts) + 1]] <- rbind(head, sub)
    hdr <- c(hdr, r + 1); r <- r + 1 + nrow(sub)
  }
  out <- do.call(rbind, parts); rownames(out) <- NULL
  list(df = out, header_rows = hdr, level_rows = setdiff(seq_len(r), hdr))
}

# -----------------------------------------------------------------------------
# Table S1. Countries
# -----------------------------------------------------------------------------
s1 <- rd("TableS1_countries.csv")
ft1 <- save_docx_table(
  s1,
  "Table S1. Countries of the included patients and classification of brucellosis endemicity",
  c(sprintf("%d patients from %d countries, ordered by number of patients; percentages are of all included patients.",
            sum(as.integer(s1[["Patients, n"]])), nrow(s1)),
    "Current status: classification used in all models (current endemicity of human brucellosis). Evidence level: A, country-specific surveillance data or official animal-health status; B, region-level statement or isolated case reports only.",
    "Classification rule: EU countries by the official brucellosis-free status of livestock; other countries by whether locally acquired human cases are reported regularly (sporadic or imported cases only = non-endemic). Reported incidence was used as supporting evidence only, because under-reporting varies widely between countries. Incidence figures are from OIE-WAHIS as compiled by Liu et al., Emerg Microbes Infect 2024, unless stated otherwise.",
    "Pappas 2006: annual incidence of human brucellosis per million population and classification from Pappas et al., Lancet Infect Dis 2006; used in the sensitivity analysis in Table S13. Status differs: the two classifications disagree."),
  file.path(OUT_TAB, "TableS1_countries.docx"),
  widths = c(1.25, 0.6, 0.6, 0.85, 0.65, 3.75, 0.95, 0.85, 0.7))
# long text reads better left-aligned
saveRDS(align(ft1, j = 6, align = "left", part = "body"),
        file.path(OUT_TAB, "TableS1_countries.rds"))

# -----------------------------------------------------------------------------
# Table S3. JBI items
# -----------------------------------------------------------------------------
s2 <- rd("TableS3_JBI_items.csv")
p2 <- stack_panels(s2, "group", c("item", "Yes", "Unclear", "No", "Not applicable"))
names(p2$df)[1] <- "JBI checklist item"
save_docx_table(
  p2$df,
  "Table S3. Reporting quality of the source reports: JBI critical appraisal checklist for case reports",
  c("Values are n (%) of patients. Each patient was rated on the eight items of the JBI critical appraisal checklist for case reports, based on the report from which the patient's data were extracted.",
    "Case reports: publications describing a single patient; case series: publications describing two or more patients."),
  file.path(OUT_TAB, "TableS3_JBI_items.docx"),
  header_rows = p2$header_rows, level_rows = p2$level_rows,
  widths = c(4.2, 1.2, 1.2, 1.2, 1.2))

# -----------------------------------------------------------------------------
# Table S6. Model-based vs CR2 standard errors
# -----------------------------------------------------------------------------
s3 <- rd("TableS6_SE_comparison.csv")
names(s3)[1] <- "Covariate"
save_docx_table(
  s3,
  "Table S6. Main model: confidence intervals with model-based and cluster-robust (CR2) standard errors",
  c(sprintf("Main log-normal accelerated failure time model, complete cases (n = %d patients from %d reports). TR, time ratio.",
            n_cc, k_cc),
    "Model-based: ordinary least-squares standard errors, which treat patients as independent. CR2: cluster-robust standard errors with Satterthwaite degrees of freedom, clustering patients by source report; used for all inference in the main analysis.",
    "SE ratio above 1: ignoring the clustering of patients within reports would understate the uncertainty of that estimate."),
  file.path(OUT_TAB, "TableS6_SE_comparison.docx"),
  widths = c(2.6, 0.8, 1.5, 1.1, 1.5, 1.1, 1.6))

# -----------------------------------------------------------------------------
# Table S5. Collinearity (the era x endemic cross-table goes into the note)
# -----------------------------------------------------------------------------
s4  <- rd("TableS5_collinearity.csv")
gap <- which(s4$Term == "")[1]
vif_tab <- s4[seq_len(gap - 1), ]
for (j in 2:3) vif_tab[[j]] <- sprintf("%.2f", as.numeric(vif_tab[[j]]))
names(vif_tab) <- c("Covariate", "VIF", "GVIF^(1/(2×Df))")
xt  <- s4[(gap + 2):nrow(s4), ]
cnt <- setNames(xt$VIF, xt$Term)
save_docx_table(
  vif_tab,
  "Table S5. Main model: collinearity between covariates",
  c(sprintf("Main model, complete cases (n = %d). VIF, variance inflation factor; GVIF, generalised VIF. Each term has one degree of freedom, so GVIF = VIF and GVIF^(1/(2×Df)) = √VIF. VIF above 5 (GVIF^(1/(2×Df)) above 2.24) would indicate problematic collinearity.", n_cc),
    sprintf("Year of publication and endemic status were the most strongly related covariates: before 2000, %s non-endemic and %s endemic patients; 2000 or later, %s non-endemic and %s endemic patients.",
            cnt["Before 2000 / Non-endemic"], cnt["Before 2000 / Endemic"],
            cnt["2000 or later / Non-endemic"], cnt["2000 or later / Endemic"])),
  file.path(OUT_TAB, "TableS5_collinearity.docx"),
  widths = c(3.6, 1.2, 1.6))

# -----------------------------------------------------------------------------
# Table S7. Influential patients
# -----------------------------------------------------------------------------
s5 <- rd("TableS7_influential_patients.csv")
dl <- as.numeric(s5$Delay)
s5$Delay <- ifelse(dl < 1, sprintf("%.2f", dl), sprintf("%.1f", dl))
s5[[4]] <- sprintf("%.2f", as.numeric(s5[[4]]))   # standardised residual
s5[[5]] <- sprintf("%.3f", as.numeric(s5[[5]]))   # Cook's distance
s5[[6]] <- sprintf("%.3f", as.numeric(s5[[6]]))   # leverage
names(s5) <- c("Patient", "Source report", "Delay, months", "Standardised residual",
               "Cook's distance", "Leverage")
save_docx_table(
  s5,
  "Table S7. Main model: influential patients. A. Patients with Cook's distance above 4/n",
  c(sprintf("Main model, complete cases (n = %d); threshold 4/n = %.3f. Patients are ordered by Cook's distance (Fig. S2C). None was removed from any analysis; panel B shows the main model refitted without them, for information only.",
            n_cc, 4 / n_cc),
    "Standardised residual: observed minus fitted log delay, divided by its standard error; positive values mean a longer delay than predicted by the model. Leverage: hat value."),
  file.path(OUT_TAB, "TableS7_influential_patients.docx"),
  widths = c(2.0, 2.0, 1.3, 1.6, 1.4, 1.2))

s5b <- rd("TableS7b_without_influential.csv")
save_docx_table(
  s5b,
  sprintf("Table S7 (continued). B. Main model with and without the %d influential patients", nrow(s5)),
  c("Values are time ratios (TR) with 95% CIs from the log-normal accelerated failure time model; CR2 cluster-robust standard errors clustered by source report.",
    "For information only: the influential patients were kept in all analyses."),
  file.path(OUT_TAB, "TableS7b_without_influential.docx"),
  widths = c(3.6, 2.2, 2.8))

# -----------------------------------------------------------------------------
# Table S8. All estimates of the secondary (logistic) models
# -----------------------------------------------------------------------------
s6 <- rd("TableS8_secondary_covariates.csv")
p6 <- stack_panels(s6, "Outcome", c("Term", "OR (95% CI)", "P"))
names(p6$df) <- c("Covariate", "Adjusted OR (95% CI)", "P")
save_docx_table(
  p6$df,
  "Table S8. Secondary outcomes: all estimates of the adjusted logistic regression models",
  c("OR, odds ratio from logistic regression. Diagnostic delay was entered as log-transformed delay and is expressed per doubling of delay (for example, from 3 to 6 months).",
    "95% CIs and P values use cluster-robust (CR2) standard errors, clustering patients by source report. Patients in whom the complication was not reported were excluded from that model (Table 3)."),
  file.path(OUT_TAB, "TableS8_secondary_covariates.docx"),
  header_rows = p6$header_rows, level_rows = p6$level_rows,
  widths = c(3.6, 2.0, 0.8))

# -----------------------------------------------------------------------------
# Collect Tables S1-S21 (sorted by table number)
# -----------------------------------------------------------------------------
TABLES <- c("TableS1_countries", "TableS3_JBI_items", "TableS6_SE_comparison",
            "TableS5_collinearity", "TableS7_influential_patients",
            "TableS8_secondary_covariates", "TableS9_report_type", "TableS10_quantile",
            "TableS11_distributions", "TableS12_jbi", "TableS13_classification_year",
            "TableS16_no_fever", "TableS15_fever_imputation", "TableS14_leave_one_out",
            "TableS17_variable_selection", "TableS18_mixed_model",
            "TableS19_report_weighted", "TableS20_period_interaction",
            "TableS2_extraction_agreement", "TableS4_prior_diagnoses",
            "TableS21_neuro_unplanned")
TABLES <- TABLES[order(as.integer(sub("^TableS(\\d+)_.*", "\\1", TABLES)))]
stopifnot(identical(as.integer(sub("^TableS(\\d+)_.*", "\\1", TABLES)), seq_along(TABLES)))
rds <- file.path(OUT_TAB, paste0(TABLES, ".rds"))
miss <- !file.exists(rds)
if (any(miss)) stop("No .rds for: ", paste(TABLES[miss], collapse = ", "),
                    "\nUpdate save_docx_table() in 00_helpers.R and rerun the script that makes each of these tables.")
csv <- sub("\\.rds$", ".csv", rds)
old <- file.exists(csv) & file.mtime(rds) < file.mtime(csv)
if (any(old)) warning("The .rds is older than the .csv (rerun the script): ",
                      paste(TABLES[old], collapse = ", "))
fts <- lapply(rds, readRDS)
for (i in seq_along(fts)) {
  cap <- paste(unlist(fts[[i]]$caption$value), collapse = " ")
  if (!grepl(sprintf("Table S%d\\.", i), cap))
    warning(TABLES[i], ": caption does not start with 'Table S", i, ".' - check the numbering in its script")
}

# -----------------------------------------------------------------------------
# Figure legends
# -----------------------------------------------------------------------------
age_kn <- round(quantile(cc$Age, c(0.05, 0.35, 0.65, 0.95), names = FALSE))
yr_kn  <- round(quantile(cc$Pub_Year, c(0.05, 0.35, 0.65, 0.95), names = FALSE))
FIGS <- list(
  list(file = "FigS2_diagnostics.png", w = 8, h = 6.5, label = "Fig. S2",
       text = sprintf(" Diagnostics of the main model (log-normal accelerated failure time model, n = %d). A, normal Q–Q plot of standardised residuals; B, standardised residuals against fitted log delay, with a loess smoother; C, Cook's distance for each patient, with the dashed line at 4/n (%.3f) and the five most influential patients labelled. Patients above the line are listed in Table S7; none was removed.",
                      n_cc, 4 / n_cc)),
  list(file = "FigS1_age_spline.png", w = 6, h = 4, label = "Fig. S1",
       text = sprintf(" Age modelled with a restricted cubic spline in the main model (4 knots at %s years: 5th, 35th, 65th and 95th percentiles). Solid line, time ratio relative to the median age (%s years) with 95%% CI from CR2 standard errors (shaded); dashed line, the linear age term used in the main model; ticks, observed ages (n = %d).",
                      paste(age_kn, collapse = ", "), round(median(cc$Age)), n_cc)),
  list(file = "FigS4_year_spline.png", w = 6, h = 4, label = "Fig. S4",
       text = sprintf(" Year of publication modelled with a restricted cubic spline (4 knots at %s: 5th, 35th, 65th and 95th percentiles), sensitivity analysis in Table S13. Solid line, time ratio relative to 2000 with 95%% CI from CR2 standard errors (shaded); dotted line, the 2000 cut-off used in the main model; ticks, observed years of publication.",
                      paste(yr_kn, collapse = ", "))),
  list(file = "FigS3_leave_one_out.png", w = 9, h = 4.5, label = "Fig. S3",
       text = " Omitting large case series (5 or more patients) in turn (Table S14). Time ratios (95% CI) for fever and for publication in 2000 or later in the main model (A), and for publication in 2000 or later in the model without fever that includes all patients (B). Each row shows the report omitted; square, no report omitted; surgical cohorts, the series that enrolled only surgically treated patients, omitted together; all large series, every report with 5 or more patients omitted. Dashed line, time ratio = 1.")
)
FIGS <- FIGS[order(vapply(FIGS, function(f) as.integer(sub("Fig. S", "", f$label, fixed = TRUE)), 0L))]   # Fig. S1-S4 in order
fig_path <- file.path(OUT_FIG, vapply(FIGS, `[[`, "", "file"))
if (any(!file.exists(fig_path))) stop("Figure not found: ",
                                      paste(fig_path[!file.exists(fig_path)], collapse = ", "))

# -----------------------------------------------------------------------------
# Build the document: A4 landscape, one table or figure per page
# -----------------------------------------------------------------------------
MAX_H <- 6.0   # inches available for a figure above its legend
f_lab <- fp_text(font.family = "Times New Roman", font.size = 10, bold = TRUE)
f_txt <- fp_text(font.family = "Times New Roman", font.size = 10)

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(left = 0.75, right = 0.75, top = 0.75, bottom = 0.75)))

for (i in seq_along(fts)) {
  # keep the footnotes together with the last table rows, so they never
  # start a new page on their own under a repeated header
  doc <- body_add_flextable(doc, paginate(fts[[i]], init = FALSE, hdr_ftr = TRUE),
                            align = "left")
  if (TABLES[i] == "TableS7_influential_patients") {   # panel B on the same page
    doc <- body_add_par(doc, "")
    doc <- body_add_flextable(doc, readRDS(file.path(OUT_TAB, "TableS7b_without_influential.rds")),
                              align = "left")
  }
  doc <- body_add_break(doc)
}
for (i in seq_along(FIGS)) {
  f <- FIGS[[i]]
  s <- min(1, MAX_H / f$h, PAGE_WIDTH / f$w)
  doc <- body_add_img(doc, src = fig_path[i], width = f$w * s, height = f$h * s)
  doc <- body_add_fpar(doc, fpar(ftext(f$label, f_lab), ftext(f$text, f_txt)))
  if (i < length(FIGS)) doc <- body_add_break(doc)
}

out_file <- file.path("outputs", "Supplementary_Material.docx")
print(doc, target = out_file)
cat("\nSaved", out_file, "with", length(fts), "tables and", length(FIGS), "figures\n")
