# =============================================================================
# 13_S7_fever_imputation.R
# Sensitivity analysis S7 (Table S15): does missing individual fever data change
# the conclusions? Fever is imputed for patients from series that reported fever
# only as a count for the whole series, and the main model is refitted in all patients.
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS15_fever_imputation.docx / .csv
#
# Allocation-based multiple imputation (m = 50):
#   in each imputation, the number of febrile patients reported by the series is
#   allocated at random among that series' patients in our data. If only part of
#   the series is in our data, the number allocated is drawn from the
#   hypergeometric distribution. The main model (CR2) is fitted to each completed
#   data set and the results are pooled with Rubin's rules (Barnard-Rubin df).
# Series-level fever counts come from the original reports (FEVER_SERIES below).
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))

M_IMP <- 50
set.seed(20261015)

# Fever counts reported for the whole series (from the original papers)
FEVER_SERIES <- data.frame(
  Cluster_ID = c("Abulizi_2017"),
  n_fever    = c(27),     # febrile patients reported in the series
  n_series   = c(32),     # patients in the series (denominator of that count)
  stringsAsFactors = FALSE)

ipd <- load_ipd()
miss <- is.na(ipd$Fever_Status)
cl_miss <- unique(ipd$Cluster_ID[miss])
not_listed <- setdiff(cl_miss, FEVER_SERIES$Cluster_ID)
if (length(not_listed)) stop("Fever missing in reports not listed in FEVER_SERIES: ",
                             paste(not_listed, collapse = ", "))
for (i in seq_len(nrow(FEVER_SERIES))) {
  n_in <- sum(miss & ipd$Cluster_ID == FEVER_SERIES$Cluster_ID[i])
  cat(sprintf("%s: %d/%d febrile in the series; %d patients with missing fever in our data\n",
              FEVER_SERIES$Cluster_ID[i], FEVER_SERIES$n_fever[i], FEVER_SERIES$n_series[i], n_in))
  if (n_in > FEVER_SERIES$n_series[i]) stop("More patients in the data than in the series")
}

impute_once <- function(d) {
  for (i in seq_len(nrow(FEVER_SERIES))) {
    idx <- which(is.na(d$Fever_Status) & d$Cluster_ID == FEVER_SERIES$Cluster_ID[i])
    n_in <- length(idx); N <- FEVER_SERIES$n_series[i]; K <- FEVER_SERIES$n_fever[i]
    k <- if (n_in == N) K else rhyper(1, K, N - K, n_in)
    v <- rep(0L, n_in); v[sample.int(n_in, k)] <- 1L
    d$Fever_Status[idx] <- v
  }
  d
}

fits <- lapply(seq_len(M_IMP), function(m) {
  d <- impute_once(ipd)
  f <- lm(main_formula(), data = d)
  t <- robust_tab(f, d$Cluster_ID)
  t <- t[match(COVARS, t$term), ]
  list(beta = setNames(t$beta, COVARS), var = setNames(t$se^2, COVARS),
       df = setNames(t$df, COVARS))
})

# Rubin's rules with the Barnard-Rubin small-sample degrees of freedom
B_mat <- sapply(fits, `[[`, "beta"); U_mat <- sapply(fits, `[[`, "var"); D_mat <- sapply(fits, `[[`, "df")
qbar <- rowMeans(B_mat); ubar <- rowMeans(U_mat); bvar <- apply(B_mat, 1, var)
tvar <- ubar + (1 + 1 / M_IMP) * bvar
lambda <- pmin(pmax((1 + 1 / M_IMP) * bvar / tvar, 1e-8), 1 - 1e-8)
nu_com <- rowMeans(D_mat)
nu_old <- (M_IMP - 1) / lambda^2
nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lambda)
nu <- nu_old * nu_obs / (nu_old + nu_obs)
se <- sqrt(tvar)
pooled <- data.frame(term = COVARS, beta = qbar, lo = qbar - qt(0.975, nu) * se,
                     hi = qbar + qt(0.975, nu) * se, p = 2 * pt(-abs(qbar / se), nu),
                     fmi = lambda)

cc <- ipd[!miss, ]
t_main <- robust_tab(lm(main_formula(), data = cc), cc$Cluster_ID)

pc <- function(t, v) { p <- t$p[t$term == v]; if (length(p)) fmt_p(p) else "" }
tab <- data.frame(
  Covariate = unname(SHORT[COVARS]),
  Main = vapply(COVARS, function(v) tr_cell(t_main, v), ""),
  Main_p = vapply(COVARS, function(v) pc(t_main, v), ""),
  MI = vapply(COVARS, function(v) tr_cell(pooled, v), ""),
  MI_p = vapply(COVARS, function(v) pc(pooled, v), ""),
  FMI = ifelse(pooled$fmi < 0.005, "<0.01", sprintf("%.2f", pooled$fmi)),
  stringsAsFactors = FALSE)
hdr <- c("Covariate",
         sprintf("Main model, complete cases\nn = %d, TR (95%% CI)", nrow(cc)), "P",
         sprintf("Fever imputed\nn = %d, TR (95%% CI)", nrow(ipd)), "P",
         "Fraction of missing\ninformation")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS15_fever_imputation.csv"), row.names = FALSE)

foot <- c(
  sprintf("Fever was imputed for %d patients from %s, using the number of febrile patients reported for the whole series (%s). In each of %d imputations these patients were allocated at random; models were pooled with Rubin's rules.",
          sum(miss), paste(sub("_(\\d{4})$", " \\1", cl_miss), collapse = ", "),
          paste(sprintf("%s: %d of %d", sub("_(\\d{4})$", " \\1", FEVER_SERIES$Cluster_ID), FEVER_SERIES$n_fever, FEVER_SERIES$n_series), collapse = "; "),
          M_IMP),
  "Values are time ratios (TR) with 95% CIs from the log-normal accelerated failure time model with CR2 cluster-robust standard errors; pooled CIs use Barnard-Rubin degrees of freedom.",
  "Fraction of missing information: share of the variance of each estimate due to the imputation (0 = no effect).")
save_docx_table(out, "Table S15. Sensitivity analysis: multiple imputation of fever", foot,
                file.path(OUT_TAB, "TableS15_fever_imputation.docx"))
print(out, row.names = FALSE)
cat("\nSaved Table S15\n")
