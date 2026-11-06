# =============================================================================
# 18_S11_report_weighted.R
# Supplementary analysis S11 (Table S19), added after the pre-specified plan:
# do large case series dominate because each patient counts equally?
# Each source report is given the same total weight (weight = 1 / number of
# its patients in the analysis), so every report counts once.
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS19_report_weighted.docx / .csv
#
# Four fits of the log-normal AFT model, CR2 cluster-robust CIs:
#   1  main model, unweighted (Table 2)
#   2  main model, equal weight per report
#   3  model without fever, all patients, unweighted (S6)
#   4  model without fever, all patients, equal weight per report
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))

ipd <- load_ipd()
COV_NF <- setdiff(COVARS, "Fever_Status")
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS)]), ]

eq_w <- function(d) 1 / ave(rep(1, nrow(d)), d$Cluster_ID, FUN = sum)
fit_set <- function(d, covars, weighted) {
  w <- if (weighted) eq_w(d) else rep(1, nrow(d))
  d$.w <- w
  robust_tab(lm(main_formula(covars), data = d, weights = .w), d$Cluster_ID)
}
tabs <- list(
  m1 = fit_set(cc, COVARS, FALSE),
  m2 = fit_set(cc, COVARS, TRUE),
  m3 = fit_set(ipd, COV_NF, FALSE),
  m4 = fit_set(ipd, COV_NF, TRUE))

pc <- function(t, v) { p <- t$p[t$term == v]; if (length(p)) fmt_p(p) else "" }
tab <- data.frame(Covariate = unname(SHORT[COVARS]), stringsAsFactors = FALSE)
for (m in names(tabs)) {
  tab[[m]] <- vapply(COVARS, function(v) tr_cell(tabs[[m]], v), "")
  tab[[paste0(m, "_p")]] <- vapply(COVARS, function(v) pc(tabs[[m]], v), "")
}
k1 <- length(unique(cc$Cluster_ID)); k3 <- length(unique(ipd$Cluster_ID))
hdr <- c("Covariate",
         sprintf("1. Main model\nn = %d (%d reports)", nrow(cc), k1), "P",
         "2. Main model,\nequal weight per report", "P",
         sprintf("3. Without fever\nn = %d (%d reports)", nrow(ipd), k3), "P",
         "4. Without fever,\nequal weight per report", "P")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS19_report_weighted.csv"), row.names = FALSE)

foot <- c(
  "Equal weight per report: each patient weighted by 1 / (number of patients from the same report in the analysis), so that every report contributes the same total weight regardless of size.",
  "Values are time ratios (TR) with 95% CIs from the log-normal accelerated failure time model; CR2 cluster-robust standard errors clustered by source report. A dash means the term is not in that model.")
save_docx_table(out, "Table S19. Additional analysis not in the analysis plan: equal weight for each source report", foot,
                file.path(OUT_TAB, "TableS19_report_weighted.docx"))
print(out, row.names = FALSE)
cat("\nSaved Table S19\n")
