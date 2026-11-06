# =============================================================================
# 11_S5_classification_year.R
# Sensitivity analysis S5 (Table S13, Fig. S4): do the endemic classification and
# the 2000 cut-off for year of publication change the conclusions?
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS13_classification_year.docx / .csv
#         outputs/figures/FigS4_year_spline.pdf / .png
# Packages: ggplot2 (plus those in 00_helpers.R)
#
# Main model (complete cases) refitted with
#   A  Pappas 2006 historical classification instead of current endemic status
#   B  year of publication as a continuous linear term (per 10 years) instead of
#      the before/after 2000 split
#   C  both changes together
#   Spline: year of publication as a restricted cubic spline (4 knots at the 5th,
#      35th, 65th and 95th percentiles) instead of the 2000 split; overall and
#      non-linearity tests; curve shown as Fig. S4
# CR2 cluster-robust CIs by source report.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages(library(ggplot2))

ipd <- load_ipd()
need <- c("Endemic_Pappas", "Pub_Year")
if (length(setdiff(need, names(ipd)))) stop("Missing: ", paste(setdiff(need, names(ipd)), collapse = ", "))
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS, need)]), ]
cc$year10 <- (cc$Pub_Year - 2000) / 10
k_rep <- length(unique(cc$Cluster_ID))
cat("Complete cases:", nrow(cc), "patients from", k_rep, "reports\n")

diff_cls <- cc$Endemic_Pappas != cc$Endemic_current
cat("Classified differently by Pappas 2006:", sum(diff_cls), "patients\n")
print(table(cc$Region[diff_cls]))
cat("Pappas 2006: endemic", sum(cc$Endemic_Pappas == 1), "| non-endemic", sum(cc$Endemic_Pappas == 0), "\n")

COV_A <- replace(COVARS, COVARS == "Endemic_current", "Endemic_Pappas")
COV_B <- replace(COVARS, COVARS == "Imaging_Era", "year10")
COV_C <- replace(COV_A, COV_A == "Imaging_Era", "year10")

fits <- list(
  main = lm(main_formula(), data = cc),
  A    = lm(main_formula(COV_A), data = cc),
  B    = lm(main_formula(COV_B), data = cc),
  C    = lm(main_formula(COV_C), data = cc))
tabs <- lapply(fits, robust_tab, cluster = cc$Cluster_ID)

# Spline for year of publication
knots <- quantile(cc$Pub_Year, c(0.05, 0.35, 0.65, 0.95), names = FALSE)
if (any(duplicated(knots))) stop("Year knots are not distinct: ", paste(knots, collapse = ", "))
cc_s <- cbind(cc, rcs_basis(cc$Pub_Year, knots, "yr_nl"))
COV_S <- c(replace(COVARS, COVARS == "Imaging_Era", "Pub_Year"), "yr_nl1", "yr_nl2")
fit_s <- lm(main_formula(COV_S), data = cc_s)
p_overall <- wald_p(fit_s, c("Pub_Year", "yr_nl1", "yr_nl2"), cc_s$Cluster_ID)
p_nonlin  <- wald_p(fit_s, c("yr_nl1", "yr_nl2"), cc_s$Cluster_ID)
cat("\nYear spline knots:", paste(round(knots), collapse = ", "),
    "| P overall", fmt_p(p_overall), "| P non-linearity", fmt_p(p_nonlin), "\n")

# Table S13
rows_terms <- c(COVARS, "year10")
row_lab <- c(SHORT[COVARS], year10 = "Publication year, per 10 years")
term_for <- function(model, v) {
  if (model %in% c("A", "C") && v == "Endemic_current") return("Endemic_Pappas")
  v
}
cell <- function(model, v) tr_cell(tabs[[model]], term_for(model, v))
pcl  <- function(model, v) {
  t <- tabs[[model]]; p <- t$p[t$term == term_for(model, v)]
  if (length(p)) fmt_p(p) else ""
}
tab <- data.frame(Covariate = unname(row_lab[rows_terms]), stringsAsFactors = FALSE)
for (m in names(fits)) {
  tab[[m]] <- vapply(rows_terms, function(v) cell(m, v), "")
  tab[[paste0(m, "_p")]] <- vapply(rows_terms, function(v) pcl(m, v), "")
}
hdr <- c("Covariate",
         "Main model\nTR (95% CI)", "P",
         "A. Pappas 2006\nclassification", "P",
         "B. Year as\ncontinuous", "P",
         "C. Both\nchanges", "P")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS13_classification_year.csv"), row.names = FALSE)

regions <- table(cc$Region[diff_cls])
foot <- c(
  sprintf("Complete cases, n = %d patients from %d reports. Values are time ratios (TR) with 95%% CIs; CR2 cluster-robust standard errors clustered by source report.", nrow(cc), k_rep),
  sprintf("A and C: endemic status by Pappas et al. (Lancet Infect Dis 2006) instead of current status; %d patients are classified differently (%s).",
          sum(diff_cls), paste(sprintf("%s %d", names(regions), as.integer(regions)), collapse = ", ")),
  "B and C: year of publication as a continuous linear term (TR per 10 years) instead of before vs after 2000; a dash means the term is not in that model.",
  sprintf("Year of publication as a restricted cubic spline (4 knots: %s): P overall = %s, P for non-linearity = %s (Fig. S4).",
          paste(round(knots), collapse = ", "), fmt_p(p_overall), fmt_p(p_nonlin)))
save_docx_table(out, "Table S13. Sensitivity analysis: endemic classification and year of publication", foot,
                file.path(OUT_TAB, "TableS13_classification_year.docx"))
print(out, row.names = FALSE)

# Fig. S4: TR across publication year, relative to 2000
yr  <- seq(min(cc$Pub_Year), max(cc$Pub_Year), length.out = 200)
ref <- 2000
L   <- cbind(yr - ref, sweep(rcs_basis(yr, knots, "yr_nl"), 2, rcs_basis(ref, knots, "yr_nl")))
tt  <- c("Pub_Year", "yr_nl1", "yr_nl2")
V   <- as.matrix(vcovCR(fit_s, cluster = cc_s$Cluster_ID, type = "CR2"))[tt, tt]
est <- drop(L %*% coef(fit_s)[tt]); se <- sqrt(rowSums((L %*% V) * L))
curve <- data.frame(year = yr, tr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))

figS3 <- ggplot(curve, aes(x = year)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = 2000, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey82") +
  geom_line(aes(y = tr), colour = "#1F1F1F", linewidth = 0.8) +
  geom_rug(data = cc, aes(x = Pub_Year), inherit.aes = FALSE, sides = "b",
           alpha = 0.4, length = unit(2, "mm")) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10), labels = c("0.5", "1", "2", "5", "10")) +   # ticks below and above 1
  # P values and knots go in the figure legend
  labs(x = "Year of publication", y = "Time ratio (vs 2000), log scale") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank())
cat("Fig. S4 legend: P overall =", fmt_p(p_overall), "; P for non-linearity =", fmt_p(p_nonlin),
    "; knots =", paste(round(knots), collapse = ", "), "\n")
print(figS3)
ggsave(file.path(OUT_FIG, "FigS4_year_spline.pdf"), figS3, width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "FigS4_year_spline.png"), figS3, width = 6, height = 4, dpi = 600)
cat("\nSaved Table S13 and Fig. S4\n")
