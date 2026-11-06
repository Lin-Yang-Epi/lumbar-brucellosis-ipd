# =============================================================================
# 12_S6_no_fever.R
# Sensitivity analysis S6 (Table S16): does leaving out the patients whose fever
# was reported only for the whole case series bias the results?
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS16_no_fever.docx / .csv
#
# Three fits of the log-normal AFT model:
#   1  main model: six covariates, complete cases (fever known)
#   2  same patients, fever removed     -> effect of dropping the fever term alone
#   3  all patients, fever removed      -> effect of adding back these patients
# Comparing 2 with 3 isolates the effect of the excluded patients.
# These patients (Abulizi_2017, 32 in the 249-patient cohort) are also described against the rest.
# CR2 cluster-robust CIs by source report.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))

ipd <- load_ipd()
COV_NF <- setdiff(COVARS, "Fever_Status")
for (v in c("Log_Delay", COV_NF)) if (anyNA(ipd[[v]])) stop(v, " has missing values")

cc  <- ipd[!is.na(ipd$Fever_Status), ]
na45 <- ipd[is.na(ipd$Fever_Status), ]
cat("Fever known:", nrow(cc), "| fever reported only for the whole series:", nrow(na45),
    "| all:", nrow(ipd), "\n")

nk <- function(d) sprintf("n = %d (%d reports)", nrow(d), length(unique(d$Cluster_ID)))
sets <- list(
  m1 = list(fit = lm(main_formula(), data = cc), d = cc),
  m2 = list(fit = lm(main_formula(COV_NF), data = cc), d = cc),
  m3 = list(fit = lm(main_formula(COV_NF), data = ipd), d = ipd))
tabs <- lapply(sets, function(s) robust_tab(s$fit, s$d$Cluster_ID))

pc <- function(m, v) { t <- tabs[[m]]; p <- t$p[t$term == v]; if (length(p)) fmt_p(p) else "" }
tab <- data.frame(Covariate = unname(SHORT[COVARS]), stringsAsFactors = FALSE)
for (m in names(sets)) {
  tab[[m]] <- vapply(COVARS, function(v) tr_cell(tabs[[m]], v), "")
  tab[[paste0(m, "_p")]] <- vapply(COVARS, function(v) pc(m, v), "")
}
hdr <- c("Covariate",
         sprintf("1. Main model\n%s", nk(cc)), "P",
         sprintf("2. Same patients,\nfever removed"), "P",
         sprintf("3. All patients,\nfever removed, %s", nk(ipd)), "P")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS16_no_fever.csv"), row.names = FALSE)

# Patients with fever reported only for the series compared with the rest
q <- function(x) { v <- quantile(x, c(0.5, 0.25, 0.75), type = 2, names = FALSE)
                   sprintf("%.1f (%.1f–%.1f)", v[1], v[2], v[3]) }
pct <- function(x) sprintf("%.1f%%", 100 * mean(x == 1))
desc <- data.frame(
  Group = c("Fever known", "Fever reported only for the series"),
  n = c(nrow(cc), nrow(na45)),
  Reports = c(length(unique(cc$Cluster_ID)), length(unique(na45$Cluster_ID))),
  Delay = c(q(cc$Delay_to_Dx), q(na45$Delay_to_Dx)),
  Before2000 = c(pct(1 - cc$Imaging_Era), pct(1 - na45$Imaging_Era)),
  Endemic = c(pct(cc$Endemic_current), pct(na45$Endemic_current)))
cat("\nPatients with fever reported only for the series (n =", nrow(na45), ") vs the rest:\n"); print(desc, row.names = FALSE)

foot <- c(
  "Values are time ratios (TR) with 95% CIs from the log-normal accelerated failure time model; CR2 cluster-robust standard errors clustered by source report. A dash means the term is not in that model.",
  "Column 2 differs from column 1 only by the fever term; column 3 differs from column 2 only by adding back the patients whose fever was reported only as a proportion for the whole case series.",
  sprintf("Patients with fever reported only for the series: n = %d from %d report%s; median delay %s months vs %s months in patients with fever known; published before 2000: %s vs %s; endemic: %s vs %s.",
          nrow(na45), desc$Reports[2], if (desc$Reports[2] == 1) "" else "s", desc$Delay[2], desc$Delay[1],
          desc$Before2000[2], desc$Before2000[1], desc$Endemic[2], desc$Endemic[1]))
save_docx_table(out, "Table S16. Sensitivity analysis: model without fever, including all patients", foot,
                file.path(OUT_TAB, "TableS16_no_fever.docx"))
print(out, row.names = FALSE)
cat("\nSaved Table S16\n")
