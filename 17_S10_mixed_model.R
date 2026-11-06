# =============================================================================
# 17_S10_mixed_model.R
# Supplementary analysis S10 (Table S18), added after the pre-specified plan:
# one-stage IPD meta-analysis with random effects, the standard IPD approach.
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS18_mixed_model.docx / .csv
# Packages: lme4, lmerTest (plus those in 00_helpers.R)
#
# Linear mixed models for log(delay) (= log-normal AFT with random effects),
# complete cases, same six covariates as the main model:
#   1  random intercept for source report
#   2  random intercepts for country and for source report (reports within countries)
# REML; fixed effects with Satterthwaite df (lmerTest). exp(coefficient) = TR.
# Variance components: share of the variance in log delay lying between reports
# (and between countries), from the null model and from the adjusted model.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

ipd <- load_ipd()
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS, "Region")]), ]
cat("Complete cases:", nrow(cc), "patients from", length(unique(cc$Cluster_ID)), "reports and",
    length(unique(cc$Region)), "countries\n")

mixed_tab <- function(fit) {
  s <- summary(fit)$coefficients
  data.frame(term = rownames(s), beta = s[, "Estimate"], se = s[, "Std. Error"],
             df = s[, "df"], p = s[, "Pr(>|t|)"],
             lo = s[, "Estimate"] - qt(0.975, s[, "df"]) * s[, "Std. Error"],
             hi = s[, "Estimate"] + qt(0.975, s[, "df"]) * s[, "Std. Error"],
             stringsAsFactors = FALSE)
}
var_share <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  v <- setNames(vc$vcov, vc$grp)
  round(100 * v / sum(v), 1)
}

RE1 <- "(1 | Cluster_ID)"
RE2 <- c("(1 | Region)", "(1 | Cluster_ID)")
fit1 <- lmer(reformulate(c(COVARS, RE1), "Log_Delay"), data = cc, REML = TRUE)
fit2 <- lmer(reformulate(c(COVARS, RE2), "Log_Delay"), data = cc, REML = TRUE)
null1 <- lmer(reformulate(c("1", RE1), "Log_Delay"), data = cc, REML = TRUE)
null2 <- lmer(reformulate(c("1", RE2), "Log_Delay"), data = cc, REML = TRUE)
for (f in list(fit1, fit2, null1, null2)) if (isSingular(f)) cat("Note: singular fit (a variance component is estimated as zero)\n")

t_main <- robust_tab(lm(main_formula(), data = cc), cc$Cluster_ID)
t1 <- mixed_tab(fit1); t2 <- mixed_tab(fit2)

s_null1 <- var_share(null1); s_fit1 <- var_share(fit1)
s_null2 <- var_share(null2); s_fit2 <- var_share(fit2)
cat("\nShare of variance in log delay (%):\n")
cat("  Report RE, null:", paste(names(s_null1), s_null1, collapse = "; "), "\n")
cat("  Report RE, adjusted:", paste(names(s_fit1), s_fit1, collapse = "; "), "\n")
cat("  Country + report RE, null:", paste(names(s_null2), s_null2, collapse = "; "), "\n")
cat("  Country + report RE, adjusted:", paste(names(s_fit2), s_fit2, collapse = "; "), "\n")

pc <- function(t, v) { p <- t$p[t$term == v]; if (length(p)) fmt_p(p) else "" }
tab <- data.frame(Covariate = unname(SHORT[COVARS]), stringsAsFactors = FALSE)
for (nm in c("main", "m1", "m2")) {
  t <- switch(nm, main = t_main, m1 = t1, m2 = t2)
  tab[[nm]] <- vapply(COVARS, function(v) tr_cell(t, v), "")
  tab[[paste0(nm, "_p")]] <- vapply(COVARS, function(v) pc(t, v), "")
}
g <- function(s, k) if (k %in% names(s)) sprintf("%.1f%%", s[[k]]) else "—"
vc_rows <- data.frame(
  Covariate = c("Variance in log delay between reports, null model",
                "Variance in log delay between reports, adjusted model",
                "Variance in log delay between countries, null model",
                "Variance in log delay between countries, adjusted model"),
  main = "—", main_p = "",
  m1 = c(g(s_null1, "Cluster_ID"), g(s_fit1, "Cluster_ID"), "—", "—"), m1_p = "",
  m2 = c(g(s_null2, "Cluster_ID"), g(s_fit2, "Cluster_ID"), g(s_null2, "Region"), g(s_fit2, "Region")),
  m2_p = "", stringsAsFactors = FALSE)
tab <- rbind(tab, vc_rows)
hdr <- c("Covariate",
         "Main model (CR2)\nTR (95% CI)", "P",
         "Random intercept:\nreport", "P",
         "Random intercepts:\ncountry and report", "P")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS18_mixed_model.csv"), row.names = FALSE)

foot <- c(
  sprintf("Complete cases, n = %d patients from %d reports in %d countries.",
          nrow(cc), length(unique(cc$Cluster_ID)), length(unique(cc$Region))),
  "One-stage IPD meta-analysis: linear mixed models of log-transformed delay (log-normal accelerated failure time model with random intercepts), REML, Satterthwaite degrees of freedom. Values are time ratios (TR) with 95% CIs.",
  "Variance rows: percentage of the total variance in log delay attributable to differences between reports (or countries), without covariates (null model) and after adjustment for the six covariates.")
if (isTRUE(all.equal(fixef(fit1), fixef(fit2), tolerance = 1e-6)))
  foot <- c(foot, "In the adjusted model the between-country variance was estimated at zero, so the model with country and report random intercepts gives the same estimates as the model with a report random intercept only.")
save_docx_table(out, "Table S18. Additional analysis not in the analysis plan: one-stage mixed-effects IPD meta-analysis", foot,
                file.path(OUT_TAB, "TableS18_mixed_model.docx"),
                first_col = 2.8)
print(out, row.names = FALSE)
cat("\nSaved Table S18\n")
