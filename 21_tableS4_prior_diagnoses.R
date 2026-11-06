# =============================================================================
# 21_tableS4_prior_diagnoses.R
# Table S4. Diagnoses and treatment before confirmation of brucellosis
# (descriptive only; companion to the 'Diagnosis before confirmation' rows of Table 1)
#
# Run from the project root, after 02_clean_data.py.
# Input : data/clean/ipd_clean.csv  (PreDx_Reported, PreDx_Category, Prior_Brucellosis)
# Output: outputs/tables/TableS4_prior_diagnoses.docx / .csv / .rds
#
# PreDx_Category holds one or more codes separated by ";" (a patient can have several
# diagnoses, so percentages can add up to more than 100%). Denominator: patients whose
# report described the period between the onset of back pain and confirmation of
# brucellosis with start of effective treatment (PreDx_Reported = 1).
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))

ipd <- load_ipd()
need <- c("PreDx_Reported", "PreDx_Category", "Prior_Brucellosis", "Endemic_current")
if (!all(need %in% names(ipd))) stop("Missing in ipd_clean.csv: ",
                                     paste(setdiff(need, names(ipd)), collapse = ", "))

CATS <- data.frame(
  code = c("D", "T", "P", "N", "R", "E", "L", "S", "V", "B1", "B0"),
  label = c("Degenerative or other non-infectious spinal disease",
            "Spinal tuberculosis",
            "Pyogenic or unspecified spinal infection",
            "Tumour, metastasis or space-occupying lesion",
            "Rheumatic or immune-mediated disease",
            "Extra-spinal disease",
            "Laboratory misidentification",
            "Symptomatic or empirical treatment only, no diagnosis named",
            "Seen by a doctor, diagnosis and treatment not reported",
            "Brucellosis diagnosed, treatment ineffective",
            "First diagnosis was brucellosis"),
  examples = c("Lumbar disc herniation, degenerative change, spondylolisthesis, spinal stenosis, osteoporotic fracture, mechanical low back pain",
               "Suspected spinal tuberculosis, or empirical antituberculous treatment",
               "Spondylodiscitis, osteomyelitis or epidural abscess without a named pathogen",
               "Vertebral metastasis, retroperitoneal or intraspinal mass",
               "Rheumatoid arthritis, psoriatic arthritis",
               "Pyelonephritis or urinary tract infection, appendicitis, amoebiasis, malaria, pneumonia, stroke or Guillain-Barré syndrome",
               "Ochrobactrum spp. or Brucella suis reported by the laboratory",
               "Analgesics, physiotherapy, self-medication with antibiotics",
               "Patient seen by a physician; the report gives no diagnosis or treatment",
               "Irregular or failed treatment, or a new complication during treatment",
               "Brucellosis diagnosed at the first visit, with no other diagnosis"),
  stringsAsFactors = FALSE)
ALT <- c("D", "T", "P", "N", "R", "E", "L")

codes <- strsplit(ifelse(is.na(ipd$PreDx_Category), "", ipd$PreDx_Category), ";")
codes <- lapply(codes, trimws)
unknown <- setdiff(unlist(codes), c(CATS$code, "U"))
if (length(unknown)) stop("Unknown PreDx_Category codes: ", paste(unknown, collapse = ", "))
is_rep <- ipd$PreDx_Reported == 1 & !is.na(ipd$PreDx_Reported)
if (any(is_rep != !vapply(codes, identical, logical(1), "U")))
  stop("PreDx_Reported does not match PreDx_Category (U = not reported)")

groups <- list(Overall = is_rep,
               Endemic = is_rep & ipd$Endemic_current == 1,
               `Non-endemic` = is_rep & ipd$Endemic_current == 0)
cell <- function(hit, g) {
  n <- sum(hit & g); d <- sum(g)
  sprintf("%d (%.1f)", n, 100 * n / d)
}
has_any <- function(keys) vapply(codes, function(x) any(x %in% keys), logical(1))

rows <- rbind(
  data.frame(Category = "Any alternative diagnosis (first seven categories)", Examples = "",
             t(vapply(groups, function(g) cell(has_any(ALT), g), "")),
             check.names = FALSE, stringsAsFactors = FALSE),
  do.call(rbind, lapply(seq_len(nrow(CATS)), function(i) {
    hit <- has_any(CATS$code[i])
    data.frame(Category = CATS$label[i], Examples = CATS$examples[i],
               t(vapply(groups, function(g) cell(hit, g), "")),
               check.names = FALSE, stringsAsFactors = FALSE)
  })))
names(rows)[3:5] <- sprintf("%s\n(n = %d)", names(groups), vapply(groups, sum, 0L))
rownames(rows) <- NULL
write.csv(rows, file.path(OUT_TAB, "TableS4_prior_diagnoses.csv"), row.names = FALSE)

n_all <- nrow(ipd); n_rep <- sum(is_rep); n_prior <- sum(ipd$Prior_Brucellosis == 1, na.rm = TRUE)
foot <- c(
  sprintf(paste("Values are n (%%). Denominator: the %d of %d patients (%.1f%%) whose report described the",
                "period between the onset of back pain and confirmation of brucellosis with start of",
                "effective treatment; for the other %d patients the report did not describe it (not reported, not taken as",
                "no other diagnosis)."), n_rep, n_all, 100 * n_rep / n_all, n_all - n_rep),
  "Diagnoses were transcribed as stated in the original reports, without inference. A patient can have more than one category, so percentages can add up to more than 100%.",
  "Brucellosis treated ineffectively: the report states irregular or failed treatment, or a new complication or lesion during treatment; persistent pain alone was not counted.",
  sprintf("Brucellosis diagnosed before the onset of back pain (%d patients) is past history and is not counted here.", n_prior))
save_docx_table(rows,
                "Table S4. Diagnoses and treatment recorded before confirmation of brucellosis",
                foot, file.path(OUT_TAB, "TableS4_prior_diagnoses.docx"),
                header_rows = integer(0), level_rows = 2:8,   # the seven categories counted as alternative diagnoses
                bold_rows = c(1, 9:nrow(rows)),   # top-level rows; 2-8 are its sub-categories
                widths = c(2.7, 4.2, 1.1, 1.1, 1.1))
print(rows, row.names = FALSE)
cat("\nSaved Table S4\n")
