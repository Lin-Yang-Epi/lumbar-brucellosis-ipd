# =============================================================================
# 22_tableS21_neuro_unplanned.R
# Unplanned analyses for the secondary outcome neurological deficit (Table S21)
#
#   Row 1  Strict definition (main definition), logistic regression as in Table 3:
#          OR per doubling of delay, CR2 cluster-robust SEs by source report,
#          restricted cubic spline (3 knots) overall and non-linearity tests
#   Row 2  Strict definition, Firth penalised logistic regression (logistf):
#          few events relative to the model degrees of freedom; profile penalised-
#          likelihood 95% CIs and penalised likelihood-ratio P values.
#          logistf has no cluster-robust option, so these CIs do not account for
#          clustering of patients within reports (stated in the table note).
#   Row 3  Broad definition (strict definition or any objective neurological sign),
#          same model as row 1
#
# Strict and broad definitions: see 02_clean_data.py (Neuro_Deficit,
# Neuro_Deficit_Broad) and Methods 2.5 / 2.8.
#
# Run from the project root, after 02_clean_data.py.
# Input : data/clean/ipd_clean.csv
# Output: outputs/tables/TableS21_neuro_unplanned.docx / .csv / .rds
# Packages: logistf (plus those in 00_helpers.R)
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages(library(logistf))

ipd <- load_ipd()
COVS2 <- c("age10", "Sex", "Endemic_current", "Imaging_Era")   # as in Table 3
need <- c("Cluster_ID", "Log_Delay", COVS2, "Neuro_Deficit", "Neuro_Deficit_Broad")
miss <- setdiff(need, names(ipd))
if (length(miss)) stop("Missing columns in ipd_clean.csv: ", paste(miss, collapse = ", "),
                       " - rerun 02_clean_data.py with the current data file")
for (v in c("Neuro_Deficit", "Neuro_Deficit_Broad")) {
  bad <- !is.na(ipd[[v]]) & !(ipd[[v]] %in% c(0, 1))
  if (any(bad)) stop(v, " must be coded 0/1")
}
# the broad definition contains the strict one and has the same missing values
stopifnot(identical(is.na(ipd$Neuro_Deficit), is.na(ipd$Neuro_Deficit_Broad)),
          all(ipd$Neuro_Deficit_Broad[ipd$Neuro_Deficit %in% 1] == 1))

fmt_or <- function(b, lo, hi, k = log(2)) sprintf("%.2f (%.2f–%.2f)", exp(k * b), exp(k * lo), exp(k * hi))

# ---- conventional logistic regression, as in Table 3 ------------------------
conv <- function(y) {
  d <- ipd[!is.na(ipd[[y]]), ]
  f_un  <- glm(reformulate("Log_Delay", y), family = binomial, data = d)
  f_lin <- glm(reformulate(c("Log_Delay", COVS2), y), family = binomial, data = d)
  kn <- quantile(d$Log_Delay, c(0.10, 0.50, 0.90), names = FALSE)
  if (any(duplicated(kn))) stop(y, ": spline knots are not distinct")
  B <- rcs_basis(d$Log_Delay, kn, "ld_nl")
  d$ld <- d$Log_Delay
  d <- cbind(d, B)
  f_spl <- glm(reformulate(c("ld", colnames(B), COVS2), y), family = binomial, data = d)
  u <- robust_tab(f_un,  d$Cluster_ID); u <- u[u$term == "Log_Delay", ]
  a <- robust_tab(f_lin, d$Cluster_ID); a <- a[a$term == "Log_Delay", ]
  list(n = nrow(d), k = length(unique(d$Cluster_ID)), ev = sum(d[[y]] == 1),
       unadj = fmt_or(u$beta, u$lo, u$hi), p_unadj = fmt_p(u$p),
       adj   = fmt_or(a$beta, a$lo, a$hi), p_adj   = fmt_p(a$p),
       spline = paste(fmt_p(wald_p(f_spl, c("ld", colnames(B)), d$Cluster_ID)), "/",
                      fmt_p(wald_p(f_spl, colnames(B), d$Cluster_ID))))
}

# ---- Firth penalised logistic regression ------------------------------------
firth <- function(y) {
  d <- ipd[!is.na(ipd[[y]]), ]
  get <- function(f) {
    i <- which(names(coef(f)) == "Log_Delay")
    list(or = fmt_or(coef(f)[i], f$ci.lower[i], f$ci.upper[i]), p = fmt_p(f$prob[i]))
  }
  fu <- get(logistf(reformulate("Log_Delay", y), data = d))
  fa <- get(logistf(reformulate(c("Log_Delay", COVS2), y), data = d))
  list(n = nrow(d), k = length(unique(d$Cluster_ID)), ev = sum(d[[y]] == 1),
       unadj = fu$or, p_unadj = fu$p, adj = fa$or, p_adj = fa$p, spline = "—")
}

res <- list(
  "Strict definition, logistic regression (Table 3)"      = conv("Neuro_Deficit"),
  "Strict definition, Firth penalised logistic regression" = firth("Neuro_Deficit"),
  "Broad definition, logistic regression"                  = conv("Neuro_Deficit_Broad"))

tab <- do.call(rbind, lapply(names(res), function(nm) {
  r <- res[[nm]]
  data.frame(Analysis = nm,
             `Patients (reports)` = sprintf("%d (%d)", r$n, r$k),
             `Events, n (%)` = sprintf("%d (%.1f%%)", r$ev, 100 * r$ev / r$n),
             `Unadjusted OR (95% CI)` = r$unadj, `P ` = r$p_unadj,
             `Adjusted OR (95% CI)` = r$adj, P = r$p_adj,
             `Spline, P overall / P non-linearity` = r$spline,
             check.names = FALSE, stringsAsFactors = FALSE)
}))
print(tab, row.names = FALSE)
cat(sprintf("\nEvents per parameter of the adjusted linear model (5 parameters besides the intercept): strict %.1f, broad %.1f\n",
            res[[1]]$ev / 5, res[[3]]$ev / 5))
write.csv(tab, file.path(OUT_TAB, "TableS21_neuro_unplanned.csv"), row.names = FALSE)

save_docx_table(
  tab,
  "Table S21. Unplanned analyses: diagnostic delay and neurological deficit with penalised regression and a broad definition",
  c("OR, odds ratio per doubling of diagnostic delay (for example, from 3 to 6 months), with delay entered as log-transformed delay. Adjusted for age, sex, current endemic status and year of publication (before 2000 vs 2000 or later), as in Table 3. Patients without a reported neurological examination at diagnosis were excluded.",
    "Strict definition (main definition, Table 3): muscle weakness or cauda equina / conus medullaris syndrome on examination at diagnosis. Broad definition: the strict definition or any objective neurological sign on examination (segmental muscle atrophy, reduced or absent tendon or superficial reflexes, upper motor neuron signs, or objective sensory loss).",
    "Logistic regression: 95% CIs and P values from cluster-robust (CR2) standard errors, clustering patients by source report; restricted cubic spline with 3 knots at the 10th, 50th and 90th percentiles of log delay. Firth penalised logistic regression: profile penalised-likelihood 95% CIs and penalised likelihood-ratio P values; these do not account for clustering of patients within reports."),
  file.path(OUT_TAB, "TableS21_neuro_unplanned.docx"),
  first_col = 2.6)
cat("Saved Table S21\n")
