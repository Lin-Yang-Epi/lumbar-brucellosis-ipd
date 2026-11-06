# =============================================================================
# 14_S8_leave_one_out.R
# Sensitivity analysis S8 (Table S14, Fig. S3): does any single large case series
# drive the results?
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS14_leave_one_out.docx / .csv
#         outputs/figures/FigS3_leave_one_out.pdf / .png
# Packages: ggplot2 (plus those in 00_helpers.R)
#
# Large series = reports contributing at least MIN_SERIES patients.
#   A  Main model (fever known): refitted leaving out each large series in turn,
#      and leaving out all of them together
#   B  Model without fever, all patients (S6 model 3): the same, so that the
#      series without individual fever data (Abulizi_2017) is also assessed
# CR2 cluster-robust CIs by source report.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages(library(ggplot2))

MIN_SERIES <- 5
# Series that enrolled only surgically treated patients (checked in the original papers);
# omitting them together is a post hoc analysis (selection on the outcome)
SURGICAL <- c("Abulizi_2017", "Li_2026")

ipd <- load_ipd()
COV_NF <- setdiff(COVARS, "Fever_Status")
size <- table(ipd$Cluster_ID)
large <- names(sort(size[size >= MIN_SERIES], decreasing = TRUE))
cat("Reports with >=", MIN_SERIES, "patients:", length(large), "\n")
print(data.frame(report = large, n = as.integer(size[large]),
                 year = ipd$Pub_Year[match(large, ipd$Cluster_ID)],
                 country = ipd$Region[match(large, ipd$Cluster_ID)],
                 fever_missing = vapply(large, function(r) sum(is.na(ipd$Fever_Status[ipd$Cluster_ID == r])), 0L)),
      row.names = FALSE)

# Report names shown as first author and year (Cluster_ID "Ariza_1985" -> "Ariza 1985")
disp <- function(x) sub("_(\\d{4})$", " \\1", x)

run_panel <- function(d, covars, panel) {
  fit_one <- function(dd) {
    cv <- varying(dd, covars)
    robust_tab(lm(main_formula(cv), data = dd), dd$Cluster_ID)
  }
  present <- large[large %in% d$Cluster_ID]
  sets <- c(list(None = d),
            setNames(lapply(present, function(r) d[d$Cluster_ID != r, ]), present),
            list(`Surgical cohorts (post hoc)` = d[!d$Cluster_ID %in% SURGICAL, ]),
            list(`All large series` = d[!d$Cluster_ID %in% present, ]))
  rows <- lapply(names(sets), function(s) {
    dd <- sets[[s]]; t <- fit_one(dd)
    lab <- if (s == "None") sprintf("None (n = %d)", nrow(dd))
           else if (s == "All large series")
             sprintf("All large series in this panel (%d reports; n = %d)", length(present), nrow(dd))
           else if (s == "Surgical cohorts (post hoc)")
             sprintf("Surgical cohorts: %s (%d patients; post hoc; n = %d)",
                     paste(disp(SURGICAL[SURGICAL %in% d$Cluster_ID]), collapse = ", "),
                     sum(d$Cluster_ID %in% SURGICAL), nrow(dd))
           else sprintf("%s (%d patients; n = %d)", disp(s), sum(d$Cluster_ID == s), nrow(dd))
    est <- lapply(covars, function(v) { r <- t[t$term == v, ]
      if (nrow(r)) c(exp(r$beta), exp(r$lo), exp(r$hi)) else c(NA, NA, NA) })
    list(label = lab, cells = vapply(covars, function(v) tr_cell(t, v), ""),
         est = setNames(est, covars), omitted = s, panel = panel)
  })
  rows
}

dA <- ipd[!is.na(ipd$Fever_Status), ]
A <- run_panel(dA, COVARS, "A")
B <- run_panel(ipd, COV_NF, "B")

mk_tab <- function(rows, covars, title) {
  body <- do.call(rbind, lapply(rows, function(r)
    data.frame(Omitted = r$label, t(r$cells), check.names = FALSE, stringsAsFactors = FALSE)))
  head <- as.data.frame(matrix("", 1, ncol(body), dimnames = list(NULL, names(body))),
                        stringsAsFactors = FALSE)
  head$Omitted <- title
  rbind(head, body)
}
tA <- mk_tab(A, COVARS, "A. Main model (patients with fever known)")
tB <- mk_tab(B, COV_NF, "B. Model without fever (all patients)")
tB$Fever_Status <- c("", rep("—", nrow(tB) - 1)); tB <- tB[, names(tA)]   # not in the model
tab <- rbind(tA, tB)
names(tab) <- c("Report omitted", unname(SHORT[COVARS]))
write.csv(tab, file.path(OUT_TAB, "TableS14_leave_one_out.csv"), row.names = FALSE)

hdr_rows <- which(tab$`Report omitted` %in% c(tA$Omitted[1], tB$Omitted[1]))
foot <- c(
  sprintf("Large series: reports contributing at least %d patients (%d reports: %s). Each was omitted in turn; the last row of each panel omits all large series included in that panel.",
          MIN_SERIES, length(large), paste(disp(large), collapse = ", ")),
  "Values are time ratios (TR) with 95% CIs from the log-normal accelerated failure time model; CR2 cluster-robust standard errors clustered by source report. A dash means the covariate did not vary after omission or is not in the model.",
  "Panel B omits the fever term so that reports without individual fever data can be included.",
  sprintf("Surgical cohorts (post hoc): series that enrolled only surgically treated patients (%s); omitted together because enrolment depended on surgery, whose indication is related to disease duration and complications, so that these series may over-represent patients with longer delays.",
          paste(disp(SURGICAL), collapse = ", ")))
save_docx_table(tab, "Table S14. Sensitivity analysis: leaving out large case series", foot,
                file.path(OUT_TAB, "TableS14_leave_one_out.docx"),
                header_rows = hdr_rows, level_rows = setdiff(seq_len(nrow(tab)), hdr_rows),
                first_col = 2.6)
print(tab, row.names = FALSE)

# Fig. S3: fever and publication year across omissions
fd <- do.call(rbind, c(
  lapply(A, function(r) data.frame(panel = "A. Main model", term = c("Fever_Status", "Imaging_Era"),
                                   label = r$label, omitted = r$omitted,
                                   tr = c(r$est$Fever_Status[1], r$est$Imaging_Era[1]),
                                   lo = c(r$est$Fever_Status[2], r$est$Imaging_Era[2]),
                                   hi = c(r$est$Fever_Status[3], r$est$Imaging_Era[3]))),
  lapply(B, function(r) data.frame(panel = "B. Without fever, all patients", term = "Imaging_Era",
                                   label = r$label, omitted = r$omitted,
                                   tr = r$est$Imaging_Era[1], lo = r$est$Imaging_Era[2],
                                   hi = r$est$Imaging_Era[3]))))
fd$term <- factor(SHORT[fd$term], levels = SHORT[c("Fever_Status", "Imaging_Era")])
fd$omitted <- factor(fd$omitted, levels = rev(unique(c("None", large, "Surgical cohorts (post hoc)", "All large series"))))
fd$is_ref <- fd$omitted == "None"

figS4 <- ggplot(fd, aes(x = tr, y = omitted)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, yend = omitted), colour = "grey35", linewidth = 0.5) +
  geom_point(aes(shape = is_ref), size = 2, colour = "#1F1F1F") +
  scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 16), guide = "none") +
  # one common x axis for all panels so that intervals can be compared directly
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2), labels = c("0.25", "0.5", "1", "2")) +
  facet_grid(. ~ panel + term) +
  scale_y_discrete(labels = disp) +
  labs(x = "Time ratio (95% CI), log scale", y = "Report omitted") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 8))
print(figS4)
ggsave(file.path(OUT_FIG, "FigS3_leave_one_out.pdf"), figS4, width = 9, height = 4.5, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "FigS3_leave_one_out.png"), figS4, width = 9, height = 4.5, dpi = 600)
cat("\nSaved Table S14 and Fig. S3\n")
