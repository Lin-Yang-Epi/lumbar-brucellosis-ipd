# =============================================================================
# 16_figS2_diagnostics.R
# Fig. S2: diagnostics of the main model (log-normal AFT, complete cases)
#   A  Normal Q-Q plot of standardised residuals (log-normal assumption)
#   B  Residuals vs fitted values (constant variance, no pattern)
#   C  Cook's distance per patient (influence); checked, nobody removed
# Also: patients above the 4/n Cook's distance threshold and, for information
# only, the main model refitted without them.
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/figures/FigS2_diagnostics.pdf / .png / .tiff
#         outputs/tables/TableS7_influential_patients.csv
# Packages: ggplot2, patchwork (plus those in 00_helpers.R)
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

ipd <- load_ipd()
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS)]), ]
fit <- lm(main_formula(), data = cc)
n <- nrow(cc)

dg <- data.frame(
  Patient_ID = cc$Patient_ID, Cluster_ID = cc$Cluster_ID, Delay = cc$Delay_to_Dx,
  fitted = fitted(fit), resid = residuals(fit), std = rstandard(fit),
  cook = cooks.distance(fit), lev = hatvalues(fit), stringsAsFactors = FALSE)
dg$idx <- seq_len(n)
thr <- 4 / n
dg$flag <- dg$cook > thr

sw <- shapiro.test(dg$std)
cat(sprintf("Shapiro-Wilk on standardised residuals: W = %.3f, P = %s (n = %d; with this n,\n  small departures give small P values - judge mainly by the Q-Q plot)\n",
            sw$statistic, fmt_p(sw$p.value), n))
cat(sprintf("Cook's distance threshold 4/n = %.4f; patients above it: %d; maximum %.3f\n",
            thr, sum(dg$flag), max(dg$cook)))

infl <- dg[dg$flag, c("Patient_ID", "Cluster_ID", "Delay", "std", "cook", "lev")]
infl <- infl[order(-infl$cook), ]
infl$std <- round(infl$std, 2); infl$cook <- round(infl$cook, 3); infl$lev <- round(infl$lev, 3)
print(infl, row.names = FALSE)
write.csv(infl, file.path(OUT_TAB, "TableS7_influential_patients.csv"), row.names = FALSE)

# Main model without the flagged patients (information only)
cc2 <- cc[!dg$flag, ]
t_all <- robust_tab(fit, cc$Cluster_ID)
t_wo  <- robust_tab(lm(main_formula(), data = cc2), cc2$Cluster_ID)
cmp <- data.frame(Covariate = unname(SHORT[COVARS]),
                  All = vapply(COVARS, function(v) tr_cell(t_all, v), ""),
                  Without_influential = vapply(COVARS, function(v) tr_cell(t_wo, v), ""))
names(cmp)[2:3] <- c(sprintf("All (n = %d)", n), sprintf("Without Cook's D > 4/n (n = %d)", nrow(cc2)))
cat("\nMain model with and without influential patients (information only):\n")
print(cmp, row.names = FALSE)
write.csv(cmp, file.path(OUT_TAB, "TableS7b_without_influential.csv"), row.names = FALSE)

# ---- Figure ------------------------------------------------------------------
th <- theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(),
                                           plot.title = element_text(face = "bold", size = 9))
qq <- qqnorm(dg$std, plot.it = FALSE)
pA <- ggplot(data.frame(x = qq$x, y = qq$y), aes(x, y)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey55") +
  geom_point(size = 1.1, alpha = 0.7, colour = "#1F1F1F") +
  labs(title = "A. Normal Q-Q plot", x = "Theoretical quantiles", y = "Standardised residuals") + th
pB <- ggplot(dg, aes(fitted, std)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 1.1, alpha = 0.7, colour = "#1F1F1F") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#2F7DBF", linewidth = 0.6) +
  labs(title = "B. Residuals vs fitted", x = "Fitted log delay", y = "Standardised residuals") + th
top <- head(dg[order(-dg$cook), ], 5)
top$hj <- ifelse(top$idx > 0.9 * max(dg$idx), 1, ifelse(top$idx < 0.1 * max(dg$idx), 0, 0.5))  # keep labels inside the panel
pC <- ggplot(dg, aes(idx, cook)) +
  geom_hline(yintercept = thr, linetype = "dashed", colour = "#B03A2E") +
  geom_segment(aes(xend = idx, yend = 0), colour = "grey50", linewidth = 0.3) +
  geom_point(aes(colour = flag), size = 1) +
  geom_text(data = top, aes(label = Patient_ID, hjust = hj), size = 2.3, vjust = -0.6) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_colour_manual(values = c(`FALSE` = "#1F1F1F", `TRUE` = "#B03A2E"), guide = "none") +
  labs(title = "C. Cook's distance", x = "Patient", y = "Cook's distance") + th

figS1 <- (pA | pB) / pC + plot_layout(heights = c(1, 0.9))
print(figS1)
ggsave(file.path(OUT_FIG, "FigS2_diagnostics.pdf"), figS1, width = 8, height = 6.5, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "FigS2_diagnostics.png"), figS1, width = 8, height = 6.5, dpi = 600)
ggsave(file.path(OUT_FIG, "FigS2_diagnostics.tiff"), figS1, width = 8, height = 6.5, dpi = 600,
       compression = "lzw")
cat("\nSaved Fig. S2 and the list of influential patients\n")
