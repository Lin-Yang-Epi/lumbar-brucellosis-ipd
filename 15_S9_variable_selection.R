# =============================================================================
# 15_S9_variable_selection.R
# Sensitivity analysis S9 (Table S17): would data-driven variable selection
# reach different conclusions from the pre-specified main model?
# For comparison only; the pre-specified main model remains the primary analysis.
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS17_variable_selection.docx / .csv
# Packages: glmnet (plus those in 00_helpers.R)
#
# Candidates: the six main-model covariates plus report type (case series vs
# case report) and JBI score. Variables that may be consequences of the delay
# (complications, number of vertebrae, department of diagnosis) are not candidates.
#   1  LASSO (log delay, standardised predictors), 10-fold cross-validation with
#      folds formed by source report; lambda.min. Repeated with 100 different fold
#      splits to give the share of runs in which each variable is selected.
#   2  Backward elimination by AIC from the full candidate model; retained
#      variables refitted with CR2 CIs (post-selection CIs are optimistic).
#      Stability: share of 500 cluster-bootstrap samples in which each variable
#      is retained.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages(library(glmnet))

N_CV_REPEATS <- 100
N_BOOT       <- 500
set.seed(20261015)

ipd <- load_ipd()
ipd$series <- as.integer(ipd$Report_Type == 2)
CAND <- c(COVARS, "series", "JBI_Score")
LAB  <- c(SHORT, series = "Case series vs case report", JBI_Score = "JBI score, per point")

cc <- ipd[complete.cases(ipd[, c("Log_Delay", CAND)]), ]
k_rep <- length(unique(cc$Cluster_ID))
cat("Complete cases:", nrow(cc), "patients from", k_rep, "reports\n")
X <- as.matrix(cc[, CAND]); y <- cc$Log_Delay

report_folds <- function(cl, k = 10) {
  u <- unique(cl); f <- sample(rep_len(seq_len(k), length(u))); f[match(cl, u)]
}

# ---- 1. LASSO --------------------------------------------------------------
cv1 <- cv.glmnet(X, y, foldid = report_folds(cc$Cluster_ID), standardize = TRUE)
b_min <- coef(cv1, s = "lambda.min")[CAND, 1]
b_1se <- coef(cv1, s = "lambda.1se")[CAND, 1]
cat(sprintf("LASSO lambda.min = %.4f; lambda.1se = %.4f\n", cv1$lambda.min, cv1$lambda.1se))
cat("Selected at lambda.min:", paste(CAND[b_min != 0], collapse = ", "), "\n")
cat("Selected at lambda.1se:", if (any(b_1se != 0)) paste(CAND[b_1se != 0], collapse = ", ") else "none", "\n")

sel <- replicate(N_CV_REPEATS, {
  cv <- cv.glmnet(X, y, foldid = report_folds(cc$Cluster_ID), standardize = TRUE)
  as.numeric(coef(cv, s = "lambda.min")[CAND, 1] != 0)
})
sel_pct <- setNames(100 * rowMeans(sel), CAND)

# ---- 2. Backward elimination by AIC ----------------------------------------
full <- lm(reformulate(CAND, "Log_Delay"), data = cc)
bw <- stats::step(full, direction = "backward", trace = 0)
kept <- attr(terms(bw), "term.labels")
cat("Retained by backward AIC:", paste(kept, collapse = ", "), "\n")
t_bw <- if (length(kept)) robust_tab(bw, cc$Cluster_ID) else data.frame(term = character(0))

cl_idx <- split(seq_len(nrow(cc)), cc$Cluster_ID)
inc <- replicate(N_BOOT, {
  idx <- unlist(cl_idx[sample.int(k_rep, k_rep, replace = TRUE)], use.names = FALSE)
  d <- cc[idx, ]
  cv <- varying(d, CAND)
  # the data are written into the model call, so step() refits on this bootstrap sample
  fit_b <- eval(bquote(lm(.(reformulate(cv, "Log_Delay")), data = .(d))))
  m <- tryCatch(stats::step(fit_b, direction = "backward", trace = 0),
                error = function(e) NULL)
  if (is.null(m)) rep(NA, length(CAND)) else CAND %in% attr(terms(m), "term.labels")
})
ok <- colSums(is.na(inc)) == 0
inc_pct <- setNames(100 * rowMeans(inc[, ok, drop = FALSE]), CAND)
cat(sprintf("Bootstrap samples usable: %d of %d\n", sum(ok), N_BOOT))

# ---- Table S17 -------------------------------------------------------------
t_main <- robust_tab(lm(main_formula(), data = cc), cc$Cluster_ID)
lasso_cell <- function(v) if (b_min[v] == 0) "Not selected" else sprintf("%.2f", exp(b_min[v]))
tab <- data.frame(
  Variable = unname(LAB[CAND]),
  Main = vapply(CAND, function(v) if (v %in% COVARS) tr_cell(t_main, v) else "Not a covariate", ""),
  Lasso = vapply(CAND, lasso_cell, ""),
  LassoPct = sprintf("%.0f%%", sel_pct),
  BW = vapply(CAND, function(v) if (v %in% kept) tr_cell(t_bw, v) else "Removed", ""),
  BootPct = sprintf("%.0f%%", inc_pct),
  stringsAsFactors = FALSE)
hdr <- c("Candidate variable",
         "Pre-specified main model\nTR (95% CI)",
         "LASSO, lambda.min\nshrunken TR",
         sprintf("LASSO: selected in\n%d CV repeats", N_CV_REPEATS),
         "Backward AIC\nTR (95% CI)",
         sprintf("Backward AIC: retained\nin %d bootstraps", N_BOOT))
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS17_variable_selection.csv"), row.names = FALSE)

foot <- c(
  sprintf("Complete cases, n = %d patients from %d reports. Candidates: the six pre-specified covariates, report type and JBI score.", nrow(cc), k_rep),
  sprintf("LASSO: linear regression of log delay with standardised predictors; 10-fold cross-validation with folds formed by source report; lambda.min = %.4f. Shrunken TR = exp(penalised coefficient). %s",
          cv1$lambda.min, if (any(b_1se != 0)) paste0("Variables selected at lambda.1se: ", paste(LAB[CAND[b_1se != 0]], collapse = ", "), ".") else "No variable was selected at lambda.1se."),
  "Backward AIC: elimination from the full candidate model; retained variables shown with CR2 cluster-robust CIs, which do not account for the selection step and are therefore too narrow. Bootstrap: cluster bootstrap by source report.",
  "Data-driven selection is shown for comparison only; the pre-specified main model (Table 2) is the primary analysis.")
save_docx_table(out, "Table S17. Sensitivity analysis: data-driven variable selection", foot,
                file.path(OUT_TAB, "TableS17_variable_selection.docx"))
print(out, row.names = FALSE)
cat("\nSaved Table S17\n")
