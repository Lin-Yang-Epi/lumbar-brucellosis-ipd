# =============================================================================
# 05_table2_fig4_main_aft.R
# Primary analysis: factors associated with diagnostic delay
# Log-normal accelerated failure time (AFT) model = linear regression of log(delay)
#
# Run from the project root (the folder that contains data/ and scripts/).
# Input : data/clean/ipd_clean.csv
# Output: outputs/tables/Table2_main_AFT.docx / .csv       unadjusted and adjusted time ratios
#         outputs/tables/TableS6_SE_comparison.csv model-based vs cluster-robust (CR2) CIs
#         outputs/tables/TableS5_collinearity.csv  VIF and GVIF^(1/(2*Df)); era x endemic table
#         outputs/figures/Fig4_forest.pdf / .png / .tiff
#         outputs/figures/FigS1_age_spline.pdf / .png       age modelled with a restricted cubic spline
#         outputs/models/main_aft.rds                      fitted models and data, reused later
# Packages: clubSandwich, ggplot2, patchwork, flextable
#
# Pre-specified covariates (no data-driven selection): age, sex, exposure history,
# current endemic status, year of publication (2000 cut-off), fever during the current illness.
# Effect measure: time ratio (TR) = exp(coefficient). TR 1.5 = 50% longer delay.
# Inference: CR2 cluster-robust standard errors, clustered by source report
# (Cluster_ID), with Satterthwaite degrees of freedom.
# =============================================================================

library(clubSandwich)
library(ggplot2)
library(patchwork)
library(flextable)

N_EXPECTED <- 249
ipd <- read.csv(file.path("data", "clean", "ipd_clean.csv"),
                na.strings = "NA", stringsAsFactors = FALSE)
stopifnot(nrow(ipd) == N_EXPECTED)

need <- c("Patient_ID", "Cluster_ID", "Delay_to_Dx", "Log_Delay", "Age", "Sex",
          "Exposure_Hx", "Endemic_current", "Imaging_Era", "Fever_Status")
miss <- setdiff(need, names(ipd))
if (length(miss)) stop("Missing columns in ipd_clean.csv: ", paste(miss, collapse = ", "))
if (any(ipd$Delay_to_Dx <= 0) || max(abs(ipd$Log_Delay - log(ipd$Delay_to_Dx))) > 1e-6)
  stop("Log_Delay is not log(Delay_to_Dx); rerun 02_clean_data.py")
for (v in c("Sex", "Exposure_Hx", "Endemic_current", "Imaging_Era", "Fever_Status")) {
  bad <- !is.na(ipd[[v]]) & !(ipd[[v]] %in% c(0, 1))
  if (any(bad)) stop(v, " must be coded 0/1")
}
for (v in setdiff(need, "Fever_Status")) {
  if (anyNA(ipd[[v]])) stop(v, " has missing values; only Fever_Status may be missing")
}

ipd$age10 <- ipd$Age / 10
COVARS <- c("age10", "Sex", "Exposure_Hx", "Endemic_current", "Imaging_Era", "Fever_Status")

# Labels: reference level first, comparison level second
spec <- data.frame(
  term  = COVARS,
  label = c("Age, per 10 years", "Sex", "Animal or raw dairy exposure",
            "Current endemic status", "Year of publication", "Fever during the current illness"),
  ref   = c(NA, "Female", "No or not mentioned", "Non-endemic", "Before 2000", "No or not mentioned"),
  cmp   = c(NA, "Male", "Yes", "Endemic", "2000 or later", "Yes"),
  short = c("Age (per 10 years)", "Male vs female", "Exposure vs no exposure",
            "Endemic vs non-endemic", "Published 2000 or later vs before", "Fever vs no fever"),
  stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------
# CR2 cluster-robust results with Satterthwaite df, one row per coefficient
robust_tab <- function(fit, cluster) {
  ct <- as.data.frame(coef_test(fit, vcov = "CR2", cluster = cluster, test = "Satterthwaite"))
  df <- ct[[grep("^df", names(ct))[1]]]
  p  <- ct[[grep("^p_", names(ct))[1]]]
  data.frame(term = names(coef(fit)), beta = ct$beta, se = ct$SE, df = df, p = p,
             lo = ct$beta - qt(0.975, df) * ct$SE,
             hi = ct$beta + qt(0.975, df) * ct$SE,
             stringsAsFactors = FALSE)
}
fmt_tr  <- function(b, lo, hi) sprintf("%.2f (%.2f–%.2f)", exp(b), exp(lo), exp(hi))
fmt_p   <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
# type = 2 matches the quantiles shown by gtsummary in Table 1
med_iqr <- function(x) {
  q <- quantile(x, c(0.5, 0.25, 0.75), type = 2, names = FALSE)
  sprintf("%.1f (%.1f–%.1f)", q[1], q[2], q[3])
}

# -----------------------------------------------------------------------------
# 2. Unadjusted models (all patients with the variable recorded)
# -----------------------------------------------------------------------------
unadj <- do.call(rbind, lapply(COVARS, function(v) {
  d <- ipd[!is.na(ipd[[v]]), ]
  fit <- lm(reformulate(v, "Log_Delay"), data = d)
  r <- robust_tab(fit, d$Cluster_ID)
  r <- r[r$term == v, ]
  r$n <- nrow(d)
  r
}))

# -----------------------------------------------------------------------------
# 3. Adjusted model (complete cases)
# -----------------------------------------------------------------------------
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS)]), ]
fit_main <- lm(reformulate(COVARS, "Log_Delay"), data = cc)
adj <- robust_tab(fit_main, cc$Cluster_ID)
adj <- adj[adj$term != "(Intercept)", ]
n_cc <- nrow(cc); k_cc <- length(unique(cc$Cluster_ID))

cat("\n==== Main model (log-normal AFT) ====\n")
cat("Complete cases:", n_cc, "patients from", k_cc, "reports",
    "(excluded:", nrow(ipd) - n_cc, "with fever reported only for the whole series)\n")
cat(sprintf("Residual SD on the log scale (log-normal scale parameter): %.3f; R-squared: %.3f\n",
            summary(fit_main)$sigma, summary(fit_main)$r.squared))
print(data.frame(term = adj$term, TR_95CI = fmt_tr(adj$beta, adj$lo, adj$hi),
                 df = round(adj$df, 1), P = fmt_p(adj$p)), row.names = FALSE)

# Model-based vs CR2 intervals (shows how much clustering matters)
mb <- summary(fit_main)$coefficients[adj$term, , drop = FALSE]
mb_lo <- mb[, "Estimate"] - qt(0.975, fit_main$df.residual) * mb[, "Std. Error"]
mb_hi <- mb[, "Estimate"] + qt(0.975, fit_main$df.residual) * mb[, "Std. Error"]
se_cmp <- data.frame(
  Term = spec$short[match(adj$term, spec$term)],
  `TR` = sprintf("%.2f", exp(adj$beta)),
  `Model-based 95% CI` = sprintf("%.2f–%.2f", exp(mb_lo), exp(mb_hi)),
  `Model-based P` = fmt_p(mb[, "Pr(>|t|)"]),
  `CR2 95% CI` = sprintf("%.2f–%.2f", exp(adj$lo), exp(adj$hi)),
  `CR2 P` = fmt_p(adj$p),
  `SE ratio (CR2 / model-based)` = sprintf("%.2f", adj$se / mb[, "Std. Error"]),
  check.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Collinearity: VIF (all terms have 1 df, so GVIF = VIF and GVIF^(1/(2*Df)) = sqrt(VIF))
# -----------------------------------------------------------------------------
X <- model.matrix(fit_main)[, -1, drop = FALSE]
vif <- diag(solve(cor(X)))
coll <- data.frame(Term = spec$short[match(colnames(X), spec$term)],
                   VIF = round(vif, 2), `GVIF^(1/(2*Df))` = round(sqrt(vif), 2),
                   check.names = FALSE)
cat("\n==== Collinearity ====\n"); print(coll, row.names = FALSE)
cat("Thresholds: VIF > 5, or GVIF^(1/(2*Df)) > 2.24, would indicate a problem\n")
if (any(vif > 5)) warning("VIF above 5: ", paste(coll$Term[vif > 5], collapse = ", "))

era_end <- table(`Year of publication` = factor(cc$Imaging_Era, 0:1, c("Before 2000", "2000 or later")),
                 `Endemic status` = factor(cc$Endemic_current, 0:1, c("Non-endemic", "Endemic")))
cat("\nYear of publication x endemic status (model sample):\n"); print(era_end)
cat(sprintf("Correlation (phi) between the two: %.2f\n",
            cor(cc$Imaging_Era, cc$Endemic_current)))

# -----------------------------------------------------------------------------
# 5. Age: restricted cubic spline (4 knots at the 5th, 35th, 65th, 95th percentiles)
# -----------------------------------------------------------------------------
rcs_basis <- function(x, knots) {           # Harrell's restricted cubic spline, non-linear part
  k <- length(knots); t <- knots
  B <- matrix(NA_real_, length(x), k - 2)
  for (j in 1:(k - 2)) {
    B[, j] <- pmax(x - t[j], 0)^3 -
      pmax(x - t[k - 1], 0)^3 * (t[k] - t[j]) / (t[k] - t[k - 1]) +
      pmax(x - t[k], 0)^3 * (t[k - 1] - t[j]) / (t[k] - t[k - 1])
  }
  colnames(B) <- paste0("age_nl", 1:(k - 2))
  B / (t[k] - t[1])^2
}
knots <- quantile(cc$Age, c(0.05, 0.35, 0.65, 0.95), names = FALSE)
if (any(duplicated(knots))) stop("Age knots are not distinct: ", paste(knots, collapse = ", "))
cc_s <- cbind(cc, rcs_basis(cc$Age, knots))
fit_spline <- lm(Log_Delay ~ Age + age_nl1 + age_nl2 + Sex + Exposure_Hx +
                   Endemic_current + Imaging_Era + Fever_Status, data = cc_s)
wt <- as.data.frame(Wald_test(fit_spline, constraints = constrain_zero(c("age_nl1", "age_nl2")),
                              vcov = "CR2", cluster = cc_s$Cluster_ID, test = "HTZ"))
p_nonlin <- wt[[grep("^p", names(wt))[1]]]
p_nonlin_mb <- anova(fit_main, fit_spline)[2, "Pr(>F)"]
cat("\n==== Age non-linearity ====\n")
cat("Knots (years):", paste(round(knots, 1), collapse = ", "), "\n")
cat(sprintf("P for non-linearity: CR2 Wald (HTZ) %s; model-based F test %s\n",
            fmt_p(p_nonlin), fmt_p(p_nonlin_mb)))
cat(if (p_nonlin < 0.05) "-> Non-linear age effect: report the spline model (Fig. S age) alongside the main model\n"
    else "-> No evidence of non-linearity: age kept linear in the main model\n")

# Spline curve: TR relative to the median age, with CR2-based 95% CI
ag  <- seq(min(cc$Age), max(cc$Age), length.out = 200)
ref <- median(cc$Age)
L   <- cbind(ag - ref, sweep(rcs_basis(ag, knots), 2, rcs_basis(ref, knots)))
tt  <- c("Age", "age_nl1", "age_nl2")
V   <- as.matrix(vcovCR(fit_spline, cluster = cc_s$Cluster_ID, type = "CR2"))[tt, tt]
est <- drop(L %*% coef(fit_spline)[tt])
se  <- sqrt(rowSums((L %*% V) * L))
age_curve <- data.frame(age = ag, tr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
                        tr_lin = exp(coef(fit_main)["age10"] * (ag - ref) / 10))

p_age <- ggplot(age_curve, aes(x = age)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
  geom_line(aes(y = tr), colour = "#1F1F1F", linewidth = 0.8) +
  geom_line(aes(y = tr_lin), colour = "#2F7DBF", linewidth = 0.7, linetype = "longdash") +
  geom_rug(data = cc, aes(x = Age), inherit.aes = FALSE, sides = "b",
           alpha = 0.4, length = unit(2, "mm")) +
  scale_y_log10(breaks = c(0.25, 0.5, 0.7, 1, 1.4, 2, 4), labels = c("0.25", "0.5", "0.7", "1", "1.4", "2", "4")) +   # log-symmetric ticks around 1
  # P for non-linearity and the line key (solid = spline, dashed = linear) go in the figure legend
  labs(x = "Age, years", y = paste0("Time ratio (vs age ", round(ref), " years), log scale")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank())
cat("Fig. S1 legend: P for non-linearity =", fmt_p(p_nonlin), "\n")

# -----------------------------------------------------------------------------
# 6. Table 2
# -----------------------------------------------------------------------------
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]; u <- unadj[unadj$term == s$term, ]; a <- adj[adj$term == s$term, ]
  if (is.na(s$ref)) {                                    # continuous
    add(Characteristic = s$label, n = as.character(nrow(ipd)),
        Delay = med_iqr(ipd$Delay_to_Dx),
        Unadj = fmt_tr(u$beta, u$lo, u$hi), Unadj_P = fmt_p(u$p),
        Adj = fmt_tr(a$beta, a$lo, a$hi), Adj_P = fmt_p(a$p), level = FALSE)
  } else {
    v <- if (s$term == "age10") "Age" else s$term
    add(Characteristic = s$label, n = "", Delay = "", Unadj = "", Unadj_P = "",
        Adj = "", Adj_P = "", level = FALSE)
    x0 <- ipd$Delay_to_Dx[ipd[[v]] %in% 0]; x1 <- ipd$Delay_to_Dx[ipd[[v]] %in% 1]
    add(Characteristic = s$ref, n = as.character(length(x0)), Delay = med_iqr(x0),
        Unadj = "1 (reference)", Unadj_P = "", Adj = "1 (reference)", Adj_P = "", level = TRUE)
    add(Characteristic = s$cmp, n = as.character(length(x1)), Delay = med_iqr(x1),
        Unadj = fmt_tr(u$beta, u$lo, u$hi), Unadj_P = fmt_p(u$p),
        Adj = fmt_tr(a$beta, a$lo, a$hi), Adj_P = fmt_p(a$p), level = TRUE)
    if (s$term == "Fever_Status") {
      xna <- ipd$Delay_to_Dx[is.na(ipd$Fever_Status)]
      add(Characteristic = "Reported only for the whole series", n = as.character(length(xna)),
          Delay = med_iqr(xna), Unadj = "Not modelled", Unadj_P = "",
          Adj = "Not modelled", Adj_P = "", level = TRUE)
    }
  }
}
tab2 <- do.call(rbind, rows)

out_tab <- file.path("outputs", "tables"); out_fig <- file.path("outputs", "figures")
out_mod <- file.path("outputs", "models")
for (d in c(out_tab, out_fig, out_mod)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

hdr <- c(Characteristic = "Characteristic", n = "Patients, n",
         Delay = "Delay, months, median (IQR)",
         Unadj = "Unadjusted TR (95% CI)", Unadj_P = "P",
         Adj = "Adjusted TR (95% CI)", Adj_P = "P")
csv_out <- tab2[, names(hdr)]; names(csv_out) <- hdr
write.csv(csv_out, file.path(out_tab, "Table2_main_AFT.csv"), row.names = FALSE)

foot <- c(
  "Diagnostic delay: from the onset of back pain (the first lumbar-related symptom of the current illness) to confirmation of brucellosis and start of effective treatment.",
  "TR, time ratio: the factor by which the expected diagnostic delay is multiplied. TR > 1 indicates a longer delay and TR < 1 a shorter delay.",
  "Log-normal accelerated failure time model (linear regression of log-transformed delay). 95% CIs and P values use cluster-robust (CR2) standard errors with Satterthwaite degrees of freedom, clustering patients by source report.",
  sprintf("Unadjusted models use every patient with the variable recorded (n = %d; fever n = %d). The adjusted model includes all six variables and the %d patients with complete data from %d reports.",
          nrow(ipd), sum(!is.na(ipd$Fever_Status)), n_cc, k_cc),
  "Fever: no or not mentioned = the report stated that the patient had no fever, or did not mention fever, coded as absent. Fever reported only as a proportion for the whole case series was treated as missing and not included in the models; see Tables S15 (multiple imputation) and S16 (model without fever).",
  "IQR, interquartile range.")

lvl <- which(tab2$level); top <- which(!tab2$level)
ft <- flextable(tab2[, names(hdr)]) |>
  set_header_labels(values = as.list(hdr)) |>
  set_caption(caption = "Table 2. Factors associated with diagnostic delay: log-normal accelerated failure time model") |>
  bold(i = top, j = 1, part = "body") |>
  bold(part = "header") |>
  padding(i = lvl, j = 1, padding.left = 14, part = "body") |>
  align(j = 2:7, align = "center", part = "all") |>
  add_footer_lines(foot) |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 9, part = "all") |>
  padding(padding.top = 1, padding.bottom = 1, part = "body") |>
  autofit()
save_as_docx(ft, path = file.path(out_tab, "Table2_main_AFT.docx"))

write.csv(se_cmp, file.path(out_tab, "TableS6_SE_comparison.csv"), row.names = FALSE)
coll_out <- rbind(coll,
                  data.frame(Term = c("", "Year of publication × endemic status (model sample)",
                                      paste0(rep(rownames(era_end), each = 2), " / ",
                                             rep(colnames(era_end), 2))),
                             VIF = c(NA, NA, as.vector(t(era_end))),
                             `GVIF^(1/(2*Df))` = NA, check.names = FALSE))
write.csv(coll_out, file.path(out_tab, "TableS5_collinearity.csv"),
          row.names = FALSE, na = "")

# -----------------------------------------------------------------------------
# 7. Fig. 4: forest plot (adjusted filled, unadjusted open)
# -----------------------------------------------------------------------------
k  <- nrow(spec)
yy <- setNames(rev(seq_len(k)), spec$term)              # first variable at the top
fp <- rbind(
  data.frame(term = adj$term, model = "Adjusted", tr = exp(adj$beta),
             lo = exp(adj$lo), hi = exp(adj$hi)),
  data.frame(term = unadj$term, model = "Unadjusted", tr = exp(unadj$beta),
             lo = exp(unadj$lo), hi = exp(unadj$hi)))
fp$y <- yy[fp$term] + ifelse(fp$model == "Adjusted", 0.14, -0.14)
fp$model <- factor(fp$model, levels = c("Adjusted", "Unadjusted"))

xr <- range(c(fp$lo, fp$hi, 0.5, 2))
xlim <- c(xr[1] / 1.15, xr[2] * 1.15)
YL <- c(0.25, k + 0.85)

p_fp <- ggplot(fp, aes(y = y, colour = model, shape = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.55) +
  geom_point(aes(x = tr), size = 2.3, stroke = 0.8, fill = "white") +
  annotate("text", x = 0.95, y = 0.4, label = "Shorter delay", hjust = 1, size = 2.7, colour = "grey30") +
  annotate("text", x = 1.05, y = 0.4, label = "Longer delay", hjust = 0, size = 2.7, colour = "grey30") +
  scale_x_log10(breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.125", "0.25", "0.5", "1", "2", "4", "8")) +
  scale_y_continuous(breaks = yy, labels = spec$short[match(names(yy), spec$term)],
                     limits = YL, expand = c(0, 0)) +
  scale_colour_manual(values = c(Adjusted = "#1F1F1F", Unadjusted = "grey50"), name = NULL) +
  scale_shape_manual(values = c(Adjusted = 15, Unadjusted = 21), name = NULL) +
  coord_cartesian(xlim = xlim, clip = "off") +
  labs(x = "Time ratio (95% CI), log scale", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(colour = "grey10", size = 8.5),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

txt <- data.frame(y = yy[adj$term], tr = fmt_tr(adj$beta, adj$lo, adj$hi), p = fmt_p(adj$p))
p_txt <- ggplot(txt, aes(y = y)) +
  geom_text(aes(x = 0, label = tr), hjust = 0, size = 2.9) +
  geom_text(aes(x = 0.78, label = p), hjust = 0, size = 2.9) +
  annotate("text", x = 0, y = k + 0.6, label = "Adjusted TR (95% CI)", hjust = 0,
           size = 2.9, fontface = "bold") +
  annotate("text", x = 0.78, y = k + 0.6, label = "P", hjust = 0, size = 2.9, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  theme_void()

fig4 <- p_fp + p_txt + plot_layout(widths = c(2.5, 1.3))
print(fig4)

w <- 8; h <- 3.6
ggsave(file.path(out_fig, "Fig4_forest.pdf"),  fig4, width = w, height = h, device = cairo_pdf)  # cairo keeps the en dash
ggsave(file.path(out_fig, "Fig4_forest.png"),  fig4, width = w, height = h, dpi = 600)
ggsave(file.path(out_fig, "Fig4_forest.tiff"), fig4, width = w, height = h, dpi = 600,
       compression = "lzw")
ggsave(file.path(out_fig, "FigS1_age_spline.pdf"), p_age, width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(out_fig, "FigS1_age_spline.png"), p_age, width = 6, height = 4, dpi = 600)

saveRDS(list(fit_main = fit_main, fit_spline = fit_spline, data_cc = cc, knots = knots,
             covariates = COVARS, spec = spec),
        file.path(out_mod, "main_aft.rds"))
cat("\nSaved Table 2, Fig. 4, Fig. S (age spline), SE comparison and collinearity tables\n")
