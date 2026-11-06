# =============================================================================
# 09_S3_distributions.R
# Sensitivity analysis S3 (Table S11): do the conclusions depend on assuming a
# log-normal distribution for the delay?
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS11_distributions.docx / .csv
# Packages: survival (plus those in 00_helpers.R)
#
# The main model refitted as parametric accelerated failure time models with
# three distributions (log-normal, Weibull, log-logistic) by maximum likelihood
# (survival::survreg). Every delay ends in a diagnosis, so there is no censoring.
# exp(coefficient) is the time ratio under each distribution. Fit is compared by
# AIC (lower = better; a difference > 2 is usually taken as meaningful).
# 95% CIs use robust (sandwich) standard errors clustered by source report; CR2
# is not available for survreg, so these are the standard cluster-robust errors.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages(library(survival))

ipd <- load_ipd()
cc <- ipd[complete.cases(ipd[, c("Delay_to_Dx", COVARS)]), ]
cc$event <- 1L
cat("Complete cases:", nrow(cc), "patients from", length(unique(cc$Cluster_ID)), "reports\n")

DISTS <- c(lognormal = "Log-normal", weibull = "Weibull", loglogistic = "Log-logistic")
f_surv <- reformulate(COVARS, "Surv(Delay_to_Dx, event)")

fits <- lapply(names(DISTS), function(dn)
  survreg(f_surv, data = cc, dist = dn, robust = TRUE, cluster = Cluster_ID))
names(fits) <- names(DISTS)

summ <- lapply(fits, function(f) {
  b <- coef(f)[COVARS]
  # f$var is the robust variance (robust = TRUE); its first rows follow coef(f)
  se <- sqrt(diag(f$var))[match(COVARS, names(coef(f)))]
  data.frame(term = COVARS, beta = b, lo = b - 1.96 * se, hi = b + 1.96 * se,
             p = 2 * pnorm(-abs(b / se)))
})
aic <- vapply(fits, AIC, numeric(1))
best <- names(which.min(aic))

cat("\nAIC:", paste(sprintf("%s %.1f", DISTS, aic), collapse = "; "), "\n")
cat("Lowest AIC:", DISTS[best], "| difference from log-normal:",
    sprintf("%.1f", aic - aic["lognormal"]), "\n")

cell <- function(dn, v) { r <- summ[[dn]]; r <- r[r$term == v, ]; fmt_tr(r$beta, r$lo, r$hi) }
pcell <- function(dn, v) fmt_p(summ[[dn]]$p[summ[[dn]]$term == v])

tab <- data.frame(Covariate = unname(SHORT[COVARS]), stringsAsFactors = FALSE)
for (dn in names(DISTS)) {
  tab[[dn]] <- vapply(COVARS, function(v) cell(dn, v), "")
  tab[[paste0(dn, "_p")]] <- vapply(COVARS, function(v) pcell(dn, v), "")
}
tab <- rbind(tab, c("AIC", unlist(lapply(names(DISTS), function(dn)
  c(sprintf("%.1f%s", aic[dn], if (dn == best) " (lowest)" else ""), "")))))

hdr <- c("Covariate", unlist(lapply(DISTS, function(x) c(paste0(x, "\nTR (95% CI)"), "P"))))
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS11_distributions.csv"), row.names = FALSE)

foot <- c(
  sprintf("Complete cases, n = %d patients from %d reports. Parametric accelerated failure time models fitted by maximum likelihood without censoring.",
          nrow(cc), length(unique(cc$Cluster_ID))),
  "TR, time ratio = exp(coefficient); TR < 1 indicates a shorter delay. 95% CIs and P values use robust standard errors clustered by source report.",
  "AIC, Akaike information criterion: lower values indicate better fit; differences > 2 are usually considered meaningful. The log-normal model gives the same coefficients as the main model, with slightly different CIs (cluster-robust instead of CR2).")
save_docx_table(out, "Table S11. Sensitivity analysis: alternative distributions for the delay", foot,
                file.path(OUT_TAB, "TableS11_distributions.docx"),
                level_rows = integer(0), first_col = 2.3)
print(out, row.names = FALSE)
cat("\nSaved Table S11\n")
