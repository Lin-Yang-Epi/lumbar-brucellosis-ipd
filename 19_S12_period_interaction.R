# =============================================================================
# 19_S12_period_interaction.R
# Supplementary analysis S12 (Table S20), added after the pre-specified plan:
# the reports span almost a century. Did the associations differ between reports
# published before 2000 and those published in 2000 or later?
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS20_period_interaction.docx / .csv
#
# Main model (complete cases) with each of the other five covariates interacting
# with publication period: one global test (CR2 Wald, HTZ, as stated in Methods;
# the likelihood-ratio test ignores clustering and is printed to the console for
# internal checking only) and one P value per covariate; plus the model fitted separately in
# each period. CR2 cluster-robust CIs by source report.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))

ipd <- load_ipd()
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS)]), ]
OTHER <- setdiff(COVARS, "Imaging_Era")
int_terms <- paste0(OTHER, ":Imaging_Era")

m_main <- lm(main_formula(), data = cc)
m_int  <- lm(main_formula(extra = int_terms), data = cc)
# R orders the variables inside an interaction name by their first appearance in
# the formula (e.g. "Imaging_Era:Fever_Status"), so look the names up in the fit
int_coef <- grep(":", names(coef(m_int)), value = TRUE)
int_name <- vapply(OTHER, function(v) int_coef[grepl(paste0("(^|:)", v, "(:|$)"), int_coef)], "")
t_int  <- robust_tab(m_int, cc$Cluster_ID)
p_glob <- wald_p(m_int, unname(int_name), cc$Cluster_ID)
lr <- 2 * (as.numeric(logLik(m_int)) - as.numeric(logLik(m_main)))
p_lr <- pchisq(lr, df = length(int_terms), lower.tail = FALSE)

strata <- lapply(0:1, function(e) {
  d <- cc[cc$Imaging_Era == e, ]
  cv <- varying(d, OTHER)
  if (length(cv) < length(OTHER))
    cat(if (e == 0) "Before 2000" else "2000 or later", ": constant covariate(s) dropped:",
        paste(setdiff(OTHER, cv), collapse = ", "), "\n")
  list(tab = robust_tab(lm(main_formula(cv), data = d), d$Cluster_ID),
       n = nrow(d), k = length(unique(d$Cluster_ID)))
})
cat(sprintf("Before 2000: %d patients (%d reports); 2000 or later: %d patients (%d reports)\n",
            strata[[1]]$n, strata[[1]]$k, strata[[2]]$n, strata[[2]]$k))
cat("Global test of the", length(int_terms), "interactions: CR2 Wald P", fmt_p(p_glob),
    "| likelihood-ratio P", fmt_p(p_lr), "\n")

pc <- function(t, v) { p <- t$p[t$term == v]; if (length(p)) fmt_p(p) else "" }
tab <- data.frame(
  Covariate = unname(SHORT[OTHER]),
  Before = vapply(OTHER, function(v) tr_cell(strata[[1]]$tab, v), ""),
  Before_p = vapply(OTHER, function(v) pc(strata[[1]]$tab, v), ""),
  After = vapply(OTHER, function(v) tr_cell(strata[[2]]$tab, v), ""),
  After_p = vapply(OTHER, function(v) pc(strata[[2]]$tab, v), ""),
  P_int = vapply(OTHER, function(v) pc(t_int, int_name[[v]]), ""),
  stringsAsFactors = FALSE)
hdr <- c("Covariate",
         sprintf("Published before 2000\nn = %d (%d reports)", strata[[1]]$n, strata[[1]]$k), "P",
         sprintf("Published 2000 or later\nn = %d (%d reports)", strata[[2]]$n, strata[[2]]$k), "P",
         "P for\ninteraction")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS20_period_interaction.csv"), row.names = FALSE)

# The dash note is added only when a covariate was constant in a period and shown as "—"
dash_note <- if (any(unlist(out) == "—"))
  " A dash means the covariate did not vary in that period and could not be estimated." else ""
foot <- c(
  paste0("Values are time ratios (TR) with 95% CIs from the log-normal accelerated failure time model fitted separately in each publication period; CR2 cluster-robust standard errors clustered by source report.",
         dash_note),
  "P: CR2 Wald test of the covariate within that publication period.",
  sprintf("P for interaction: CR2 Wald test of each covariate × publication period interaction in a model with all %d interactions; global test of all %d interactions (CR2 Wald) P = %s.",
          length(int_terms), length(int_terms), fmt_p(p_glob)))
save_docx_table(out, "Table S20. Additional analysis not in the analysis plan: associations by publication period", foot,
                file.path(OUT_TAB, "TableS20_period_interaction.docx"))
print(out, row.names = FALSE)
cat("\nSaved Table S20\n")
