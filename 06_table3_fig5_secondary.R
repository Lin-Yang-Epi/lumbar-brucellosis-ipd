# =============================================================================
# 06_table3_fig5_secondary.R
# Secondary analysis: is a longer diagnostic delay associated with
# epidural abscess and neurological deficit at diagnosis?
#
# Run from the project root (the folder that contains data/ and scripts/).
# Input : data/clean/ipd_clean.csv
# Output: outputs/tables/Table3_secondary.docx / .csv
#         outputs/tables/TableS8_secondary_covariates.csv   all ORs of the adjusted models
#         outputs/figures/Fig5_delay_complications.pdf / .png / .tiff
#         outputs/models/secondary_logistic.rds
# Packages: clubSandwich, ggplot2, flextable
#
# Pre-specified: logistic regression; exposure = log(delay); adjusted for age,
# sex, current endemic status and year of publication (2000 cut-off).
#   - Linear term: OR per doubling of delay
#   - Restricted cubic spline, 3 knots (10th/50th/90th percentiles of log delay):
#     overall test and test of non-linearity
# Inference: CR2 cluster-robust standard errors by source report.
# Patients with the complication not reported are excluded from that model.
# Associations only: a complication may also appear early with virulent
# infection or impaired host defence, so no causal reading.
# =============================================================================

library(clubSandwich)
library(ggplot2)
library(flextable)

N_EXPECTED <- 249
ipd <- read.csv(file.path("data", "clean", "ipd_clean.csv"),
                na.strings = "NA", stringsAsFactors = FALSE)
stopifnot(nrow(ipd) == N_EXPECTED)

OUTCOMES <- c(Epidural_Abscess = "Epidural abscess", Neuro_Deficit = "Neurological deficit")
COVS     <- c("age10", "Sex", "Endemic_current", "Imaging_Era")
COV_LAB  <- c(age10 = "Age, per 10 years", Sex = "Male vs female",
              Endemic_current = "Endemic vs non-endemic",
              Imaging_Era = "Published 2000 or later vs before")

need <- c("Cluster_ID", "Delay_to_Dx", "Log_Delay", "Age", "Sex", "Endemic_current",
          "Imaging_Era", names(OUTCOMES))
miss <- setdiff(need, names(ipd))
if (length(miss)) stop("Missing columns in ipd_clean.csv: ", paste(miss, collapse = ", "))
for (v in c("Sex", "Endemic_current", "Imaging_Era", names(OUTCOMES))) {
  bad <- !is.na(ipd[[v]]) & !(ipd[[v]] %in% c(0, 1))
  if (any(bad)) stop(v, " must be coded 0/1")
}
for (v in c("Cluster_ID", "Delay_to_Dx", "Log_Delay", "Age", "Sex", "Endemic_current", "Imaging_Era")) {
  if (anyNA(ipd[[v]])) stop(v, " has missing values")
}
ipd$age10 <- ipd$Age / 10

# Delay groups for the observed event rates
DELAY_BREAKS <- c(-Inf, 1, 3, 6, Inf)
DELAY_LABELS <- c("≤1 month", ">1 to ≤3 months", ">3 to ≤6 months", ">6 months")
ipd$delay_grp <- cut(ipd$Delay_to_Dx, DELAY_BREAKS, labels = DELAY_LABELS, right = TRUE)

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------
robust_tab <- function(fit, cluster) {
  ct <- as.data.frame(coef_test(fit, vcov = "CR2", cluster = cluster, test = "Satterthwaite"))
  df <- ct[[grep("^df", names(ct))[1]]]
  p  <- ct[[grep("^p_", names(ct))[1]]]
  data.frame(term = names(coef(fit)), beta = ct$beta, se = ct$SE, df = df, p = p,
             lo = ct$beta - qt(0.975, df) * ct$SE,
             hi = ct$beta + qt(0.975, df) * ct$SE, stringsAsFactors = FALSE)
}
wald_p <- function(fit, terms, cluster) {
  wt <- as.data.frame(Wald_test(fit, constraints = constrain_zero(terms),
                                vcov = "CR2", cluster = cluster, test = "HTZ"))
  wt[[grep("^p", names(wt))[1]]]
}
rcs_basis <- function(x, knots, prefix) {   # Harrell's restricted cubic spline, non-linear part
  k <- length(knots); t <- knots
  B <- matrix(NA_real_, length(x), k - 2)
  for (j in 1:(k - 2)) {
    B[, j] <- pmax(x - t[j], 0)^3 -
      pmax(x - t[k - 1], 0)^3 * (t[k] - t[j]) / (t[k] - t[k - 1]) +
      pmax(x - t[k], 0)^3 * (t[k - 1] - t[j]) / (t[k] - t[k - 1])
  }
  colnames(B) <- paste0(prefix, 1:(k - 2))
  B / (t[k] - t[1])^2
}
# OR for a change of k units in the coefficient's variable (k = log(2): per doubling of delay)
fmt_or <- function(b, lo, hi, k = 1) sprintf("%.2f (%.2f–%.2f)", exp(k * b), exp(k * lo), exp(k * hi))
fmt_p  <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# -----------------------------------------------------------------------------
# 2. Models, one outcome at a time
# -----------------------------------------------------------------------------
rows <- list(); curves <- list(); observed <- list(); cov_out <- list(); fits <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

for (y in names(OUTCOMES)) {
  lab <- OUTCOMES[[y]]
  d <- ipd[!is.na(ipd[[y]]), ]
  n <- nrow(d); ev <- sum(d[[y]] == 1); k_rep <- length(unique(d$Cluster_ID))

  f_un  <- glm(reformulate("Log_Delay", y), family = binomial, data = d)
  f_lin <- glm(reformulate(c("Log_Delay", COVS), y), family = binomial, data = d)

  kn <- quantile(d$Log_Delay, c(0.10, 0.50, 0.90), names = FALSE)
  if (any(duplicated(kn))) stop(lab, ": spline knots are not distinct: ", paste(round(exp(kn), 2), collapse = ", "))
  B <- rcs_basis(d$Log_Delay, kn, "ld_nl")
  d$ld <- d$Log_Delay
  d <- cbind(d, B)
  f_spl <- glm(reformulate(c("ld", colnames(B), COVS), y), family = binomial, data = d)

  u <- robust_tab(f_un, d$Cluster_ID);  u <- u[u$term == "Log_Delay", ]
  a_all <- robust_tab(f_lin, d$Cluster_ID); a <- a_all[a_all$term == "Log_Delay", ]
  p_overall <- wald_p(f_spl, c("ld", colnames(B)), d$Cluster_ID)
  p_nonlin  <- wald_p(f_spl, colnames(B), d$Cluster_ID)

  cat("\n====", lab, "====\n")
  cat(sprintf("Patients %d (from %d reports); events %d (%.1f%%); not reported %d\n",
              n, k_rep, ev, 100 * ev / n, sum(is.na(ipd[[y]]))))
  cat(sprintf("Events per model degree of freedom (spline model, 6 df): %.1f\n", ev / 6))
  cat("Spline knots (months):", paste(round(exp(kn), 2), collapse = ", "), "\n")
  cat("OR per doubling of delay: unadjusted", fmt_or(u$beta, u$lo, u$hi, log(2)), "P", fmt_p(u$p),
      "| adjusted", fmt_or(a$beta, a$lo, a$hi, log(2)), "P", fmt_p(a$p), "\n")
  cat("Spline (adjusted): P overall", fmt_p(p_overall), "| P non-linearity", fmt_p(p_nonlin), "\n")

  # Table rows
  add(Characteristic = sprintf("%s (n = %d; %d events)", lab, n, ev), Events = "",
      Unadj = "", Unadj_P = "", Adj = "", Adj_P = "", level = FALSE)
  tb <- table(d$delay_grp, d[[y]])
  for (g in DELAY_LABELS) {
    nn <- sum(tb[g, ]); e <- if ("1" %in% colnames(tb)) tb[g, "1"] else 0
    add(Characteristic = paste("Delay", g), Events = sprintf("%d/%d (%.1f%%)", e, nn, 100 * e / nn),
        Unadj = "", Unadj_P = "", Adj = "", Adj_P = "", level = TRUE)
  }
  add(Characteristic = "Per doubling of diagnostic delay", Events = "",
      Unadj = fmt_or(u$beta, u$lo, u$hi, log(2)), Unadj_P = fmt_p(u$p),
      Adj = fmt_or(a$beta, a$lo, a$hi, log(2)), Adj_P = fmt_p(a$p), level = TRUE)
  add(Characteristic = "Restricted cubic spline, P overall / P non-linearity", Events = "",
      Unadj = "", Unadj_P = "", Adj = paste(fmt_p(p_overall), "/", fmt_p(p_nonlin)), Adj_P = "",
      level = TRUE)

  # All ORs of the adjusted linear model (supplement)
  cv <- a_all[a_all$term != "(Intercept)", ]
  cov_out[[y]] <- data.frame(
    Outcome = lab,
    Term = ifelse(cv$term == "Log_Delay", "Diagnostic delay, per doubling", COV_LAB[cv$term]),
    `OR (95% CI)` = ifelse(cv$term == "Log_Delay",
                           fmt_or(cv$beta, cv$lo, cv$hi, log(2)), fmt_or(cv$beta, cv$lo, cv$hi)),
    P = fmt_p(cv$p), check.names = FALSE)

  # Curve: average predicted probability over the sample at each delay
  # (marginal standardisation), 95% CI by the delta method with the CR2 covariance
  V <- as.matrix(vcovCR(f_spl, cluster = d$Cluster_ID, type = "CR2"))
  b <- coef(f_spl); V <- V[names(b), names(b)]
  X0 <- model.matrix(f_spl)
  lgrid <- seq(quantile(d$Log_Delay, 0.025), quantile(d$Log_Delay, 0.975), length.out = 150)
  cr <- t(vapply(lgrid, function(l) {
    X <- X0
    X[, "ld"] <- l
    X[, colnames(B)] <- rcs_basis(rep(l, nrow(X)), kn, "ld_nl")
    p <- plogis(drop(X %*% b)); m <- mean(p)
    g <- colMeans(p * (1 - p) * X)
    se_logit <- sqrt(drop(t(g) %*% V %*% g)) / (m * (1 - m))
    c(m, plogis(qlogis(m) - 1.96 * se_logit), plogis(qlogis(m) + 1.96 * se_logit))
  }, numeric(3)))
  facet <- sprintf("%s\nP overall = %s; P non-linearity = %s", lab, fmt_p(p_overall), fmt_p(p_nonlin))
  curves[[y]] <- data.frame(outcome = facet, delay = exp(lgrid),
                            prob = cr[, 1], lo = cr[, 2], hi = cr[, 3])
  ob <- data.frame(delay_grp = DELAY_LABELS,
                   event = as.vector(tapply(d[[y]], d$delay_grp, mean)),
                   n     = as.vector(table(d$delay_grp)),
                   delay = as.vector(tapply(d$Delay_to_Dx, d$delay_grp, median)),
                   outcome = facet, stringsAsFactors = FALSE)
  ob <- ob[ob$n > 0, ]
  observed[[y]] <- ob

  fits[[y]] <- list(unadjusted = f_un, adjusted_linear = f_lin, adjusted_spline = f_spl,
                    knots = kn, data = d)
}

# -----------------------------------------------------------------------------
# 3. Table 3
# -----------------------------------------------------------------------------
tab3 <- do.call(rbind, rows)
out_tab <- file.path("outputs", "tables"); out_fig <- file.path("outputs", "figures")
out_mod <- file.path("outputs", "models")
for (dd in c(out_tab, out_fig, out_mod)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

hdr <- c(Characteristic = "Outcome and diagnostic delay", Events = "Patients with the complication, n/N (%)",
         Unadj = "Unadjusted OR (95% CI)", Unadj_P = "P",
         Adj = "Adjusted OR (95% CI)", Adj_P = "P")
csv_out <- tab3[, names(hdr)]; names(csv_out) <- hdr
write.csv(csv_out, file.path(out_tab, "Table3_secondary.csv"), row.names = FALSE)
write.csv(do.call(rbind, cov_out), file.path(out_tab, "TableS8_secondary_covariates.csv"),
          row.names = FALSE)

foot <- c(
  "Diagnostic delay: from the onset of back pain (the first lumbar-related symptom of the current illness) to confirmation of brucellosis and start of effective treatment.",
  "OR, odds ratio per doubling of diagnostic delay (for example, from 3 to 6 months), from logistic regression with log-transformed delay.",
  "Adjusted for age, sex, current endemic status and year of publication (before 2000 vs 2000 or later). 95% CIs and P values use cluster-robust (CR2) standard errors, clustering patients by source report.",
  "Restricted cubic spline: delay modelled with 3 knots at the 10th, 50th and 90th percentiles of log delay in the adjusted model; P overall tests any association with delay, P non-linearity tests departure from the linear term.",
  "Patients in whom the complication was not reported (no imaging, operative or neurological examination findings described for the patient, or findings given only for the whole case series) were excluded from that analysis. Associations are not causal: a complication may also be present early in virulent infection or in patients with impaired host defence.")

lvl <- which(tab3$level); top <- which(!tab3$level)
ft <- flextable(tab3[, names(hdr)]) |>
  set_header_labels(values = as.list(hdr)) |>
  set_caption(caption = "Table 3. Diagnostic delay and complications at diagnosis") |>
  bold(i = top, j = 1, part = "body") |>
  bold(part = "header") |>
  padding(i = lvl, j = 1, padding.left = 14, part = "body") |>
  align(j = 2:6, align = "center", part = "all") |>
  add_footer_lines(foot) |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 9, part = "all") |>
  padding(padding.top = 1, padding.bottom = 1, part = "body") |>
  autofit()
save_as_docx(ft, path = file.path(out_tab, "Table3_secondary.docx"))

# -----------------------------------------------------------------------------
# 4. Fig. 5: adjusted probability of each complication across the delay range
# -----------------------------------------------------------------------------
cd <- do.call(rbind, curves); od <- do.call(rbind, observed)
lev <- unique(cd$outcome); cd$outcome <- factor(cd$outcome, lev); od$outcome <- factor(od$outcome, lev)
xb <- c(0.25, 0.5, 1, 2, 3, 6, 12, 24, 60)

fig5 <- ggplot(cd, aes(x = delay)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey82") +
  geom_line(aes(y = prob), colour = "#1F1F1F", linewidth = 0.8) +
  geom_point(data = od, aes(y = event, size = n), shape = 21, fill = "#2F7DBF",
             colour = "white", stroke = 0.4, alpha = 0.9) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_x_log10(breaks = xb, labels = format(xb, drop0trailing = TRUE, trim = TRUE)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
  scale_size_area(max_size = 6, name = "Patients in delay group") +
  labs(x = "Diagnostic delay, months (log scale)", y = "Probability of the complication") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5, lineheight = 1.1),
        legend.position = "bottom",
        panel.spacing.x = unit(8, "mm"),
        plot.margin = margin(5, 20, 5, 5))   # room for the last x label
print(fig5)

w <- 8; h <- 4
ggsave(file.path(out_fig, "Fig5_delay_complications.pdf"), fig5, width = w, height = h, device = cairo_pdf)
ggsave(file.path(out_fig, "Fig5_delay_complications.png"), fig5, width = w, height = h, dpi = 600)
ggsave(file.path(out_fig, "Fig5_delay_complications.tiff"), fig5, width = w, height = h, dpi = 600,
       compression = "lzw")

saveRDS(fits, file.path(out_mod, "secondary_logistic.rds"))
cat("\nSaved Table 3, Fig. 5 and the covariate table\n")
cat("Curve: average predicted probability (adjusted spline model); dots: observed proportion",
    "in each delay group, placed at the group's median delay\n")
