# =============================================================================
# 08_S2_quantile.R
# Sensitivity analysis S2 (Table S10): are the results driven by extreme delays?
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS10_quantile.docx / .csv
# Packages: quantreg (plus those in 00_helpers.R)
#
# Quantile regression of log(delay) on the six main-model covariates, complete cases.
# Because quantiles are preserved by the log transform, exp(coefficient) is the
# ratio of the delay quantile between groups (TR for the median, etc.).
#   - Median (50th percentile): the pre-specified S2 analysis
#   - 25th and 75th percentiles: shorter and longer delays, for context
# 95% CIs: cluster bootstrap by source report (resampling whole reports), 2,000
# replicates, percentile intervals and percentile bootstrap P values.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages(library(quantreg))

B_REPS <- 2000
TAUS   <- c(0.50, 0.25, 0.75)
set.seed(20261015)

ipd <- load_ipd()
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS)]), ]
cat("Complete cases:", nrow(cc), "patients from", length(unique(cc$Cluster_ID)), "reports\n")

t_main <- robust_tab(lm(main_formula(), data = cc), cc$Cluster_ID)

fit_rq <- function(d, tau) {
  f <- suppressWarnings(rq(main_formula(), tau = tau, data = d, method = "br"))
  coef(f)[COVARS]
}

cl_idx <- split(seq_len(nrow(cc)), cc$Cluster_ID)
k <- length(cl_idx)

res <- lapply(TAUS, function(tau) {
  est <- fit_rq(cc, tau)
  boot <- t(replicate(B_REPS, {
    idx <- unlist(cl_idx[sample.int(k, k, replace = TRUE)], use.names = FALSE)
    tryCatch(fit_rq(cc[idx, ], tau), error = function(e) rep(NA_real_, length(COVARS)))
  }))
  ok <- complete.cases(boot)
  boot <- boot[ok, , drop = FALSE]
  ci <- apply(boot, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
  # Percentile bootstrap P (consistent with the percentile CI): twice the smaller
  # share of replicates on either side of zero, at least 1/number of replicates
  p  <- apply(boot, 2, function(b) max(2 * min(mean(b <= 0), mean(b >= 0)), 1 / length(b)))
  p  <- pmin(p, 1)
  cat(sprintf("tau = %.2f: %d of %d bootstrap replicates usable\n", tau, sum(ok), B_REPS))
  data.frame(term = COVARS, beta = est, lo = ci[1, ], hi = ci[2, ], p = p, tau = tau)
})
names(res) <- paste0("q", TAUS * 100)

cell <- function(r, v) { x <- r[r$term == v, ]; fmt_tr(x$beta, x$lo, x$hi) }
pcell <- function(r, v) fmt_p(r$p[r$term == v])
tab <- data.frame(
  Covariate = unname(SHORT[COVARS]),
  Main = vapply(COVARS, function(v) tr_cell(t_main, v), ""),
  Q50 = vapply(COVARS, function(v) cell(res$q50, v), ""),
  P50 = vapply(COVARS, function(v) pcell(res$q50, v), ""),
  Q25 = vapply(COVARS, function(v) cell(res$q25, v), ""),
  Q75 = vapply(COVARS, function(v) cell(res$q75, v), ""),
  stringsAsFactors = FALSE)
hdr <- c(Covariate = "Covariate",
         Main = "Main model\nTR (95% CI)",
         Q50 = "Median regression\nratio (95% CI)", P50 = "P",
         Q25 = "25th percentile\nratio (95% CI)",
         Q75 = "75th percentile\nratio (95% CI)")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS10_quantile.csv"), row.names = FALSE)

foot <- c(
  sprintf("Complete cases, n = %d patients from %d reports.", nrow(cc), k),
  "Quantile regression of log-transformed delay; exp(coefficient) is the ratio of the delay quantile (for example, the median delay) between groups, adjusted for the other covariates. Ratio < 1 indicates a shorter delay.",
  sprintf("95%% CIs are percentile intervals from a cluster bootstrap that resamples whole source reports (%s replicates); P values are percentile bootstrap P values, consistent with the CIs.",
          format(B_REPS, big.mark = ",")),
  "The main model column repeats Table 2 (log-normal accelerated failure time model, CR2 cluster-robust CIs).")
save_docx_table(out, "Table S10. Sensitivity analysis: quantile regression", foot,
                file.path(OUT_TAB, "TableS10_quantile.docx"))
print(out, row.names = FALSE)
cat("\nSaved Table S10\n")
