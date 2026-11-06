# =============================================================================
# 10_S4_jbi.R
# Sensitivity analysis S4 (Table S12): does the reporting quality of the source
# case reports change the conclusions?
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS12_jbi.docx / .csv
#
# The main model (complete cases) refitted
#   - in patients whose report scored 8/8 on the JBI case report checklist
#     (pre-specified S4)
#   - in patients scoring 7 or 8 (less strict cut-off, for context)
# CR2 cluster-robust CIs by source report. A covariate that does not vary in a
# subset is dropped from that model and shown as a dash.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))

ipd <- load_ipd()
if (!"JBI_Score" %in% names(ipd)) stop("JBI_Score missing from ipd_clean.csv")
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS, "JBI_Score")]), ]
cat("Complete cases:", nrow(cc), "| JBI score distribution:\n"); print(table(cc$JBI_Score))

fit_subset <- function(d, label) {
  cv <- varying(d, COVARS)
  if (length(cv) < length(COVARS))
    cat(label, ": constant covariate(s) dropped:", paste(setdiff(COVARS, cv), collapse = ", "), "\n")
  t <- robust_tab(lm(main_formula(cv), data = d), d$Cluster_ID)
  list(tab = t, n = nrow(d), k = length(unique(d$Cluster_ID)))
}

sets <- list(
  main = fit_subset(cc, "All complete cases"),
  ge7  = fit_subset(cc[cc$JBI_Score >= 7, ], "JBI 7-8"),
  eq8  = fit_subset(cc[cc$JBI_Score == 8, ], "JBI 8")
)
for (s in names(sets)) cat(sprintf("%s: %d patients from %d reports\n", s, sets[[s]]$n, sets[[s]]$k))

pc <- function(s, v) { t <- sets[[s]]$tab; p <- t$p[t$term == v]; if (length(p)) fmt_p(p) else "" }
tab <- data.frame(Covariate = unname(SHORT[COVARS]), stringsAsFactors = FALSE)
for (s in names(sets)) {
  tab[[s]] <- vapply(COVARS, function(v) tr_cell(sets[[s]]$tab, v), "")
  tab[[paste0(s, "_p")]] <- vapply(COVARS, function(v) pc(s, v), "")
}
lab <- c(main = "All patients (main model)", ge7 = "JBI score 7–8", eq8 = "JBI score 8")
hdr <- c("Covariate", unlist(lapply(names(sets), function(s)
  c(sprintf("%s\nn = %d, TR (95%% CI)", lab[s], sets[[s]]$n), "P"))))
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS12_jbi.csv"), row.names = FALSE)

foot <- c(
  "JBI, Joanna Briggs Institute critical appraisal checklist for case reports (8 items; score = number of items rated Yes).",
  sprintf("Patients (reports): all %d (%d); JBI 7–8 %d (%d); JBI 8 %d (%d). Complete cases for the six covariates.",
          sets$main$n, sets$main$k, sets$ge7$n, sets$ge7$k, sets$eq8$n, sets$eq8$k),
  "TR, time ratio from the log-normal accelerated failure time model; 95% CIs and P values use CR2 cluster-robust standard errors clustered by source report. A dash means the covariate did not vary in that subset and was not estimated.")
save_docx_table(out, "Table S12. Sensitivity analysis: restriction by reporting quality", foot,
                file.path(OUT_TAB, "TableS12_jbi.docx"))
print(out, row.names = FALSE)
cat("\nSaved Table S12\n")
