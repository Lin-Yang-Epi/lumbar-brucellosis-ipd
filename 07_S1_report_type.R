# =============================================================================
# 07_S1_report_type.R
# Sensitivity analysis S1 (Table S9): do patients from single case reports differ
# from patients in case series? Single case reports are more likely to be
# selected for being unusual.
#
# Run from the project root. Input: data/clean/ipd_clean.csv
# Output: outputs/tables/TableS9_report_type.docx / .csv
#
# For each definition of report type
#   primary    : case report = 1 patient in the publication; case series = 2 or more
#   robustness : case report = 1 or 2 patients; case series = 3 or more
# the main model (complete cases) is refitted
#   (a) with report type added as a covariate
#   (b) separately in case reports and in case series
#   (c) with report type x covariate interactions: one global test (6 df) and
#       one P value per covariate
# Tests: CR2 Wald (HTZ), as stated in Methods; the likelihood-ratio test (ignores
# clustering) is printed to the console for internal checking only, not reported.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
ipd <- load_ipd()

need <- c("Report_Type", "Report_Type_alt")
if (length(setdiff(need, names(ipd)))) stop("Report_Type / Report_Type_alt missing; rerun 02_clean_data.py")
cc <- ipd[complete.cases(ipd[, c("Log_Delay", COVARS, need)]), ]
cat("Complete cases:", nrow(cc), "patients from", length(unique(cc$Cluster_ID)), "reports\n")

fit_main <- lm(main_formula(), data = cc)
t_main <- robust_tab(fit_main, cc$Cluster_ID)

run_definition <- function(var, title) {
  d <- cc
  d$series <- as.integer(d[[var]] == 2)
  n_cr <- sum(d$series == 0); n_cs <- sum(d$series == 1)
  k_cr <- length(unique(d$Cluster_ID[d$series == 0])); k_cs <- length(unique(d$Cluster_ID[d$series == 1]))

  # (a) report type as a covariate
  m_add <- lm(main_formula(extra = "series"), data = d)
  t_add <- robust_tab(m_add, d$Cluster_ID)

  # (b) stratified
  strat <- lapply(0:1, function(s) {
    ds <- d[d$series == s, ]
    cv <- varying(ds, COVARS)
    if (length(cv) < length(COVARS))
      cat("  ", title, if (s == 0) "case reports" else "case series",
          ": constant covariate(s) dropped:", paste(setdiff(COVARS, cv), collapse = ", "), "\n")
    robust_tab(lm(main_formula(cv), data = ds), ds$Cluster_ID)
  })

  # (c) interactions
  m_int <- lm(reformulate(c(COVARS, "series", paste0(COVARS, ":series")), "Log_Delay"), data = d)
  int_terms <- paste0(COVARS, ":series")
  t_int <- robust_tab(m_int, d$Cluster_ID)
  p_global_cr2 <- wald_p(m_int, int_terms, d$Cluster_ID)
  lr <- 2 * (as.numeric(logLik(m_int)) - as.numeric(logLik(m_add)))
  p_global_lr <- pchisq(lr, df = length(int_terms), lower.tail = FALSE)

  cat("\n====", title, "====\n")
  cat(sprintf("Case reports %d patients (%d reports); case series %d patients (%d reports)\n",
              n_cr, k_cr, n_cs, k_cs))
  cat("Case series vs case report, TR:", tr_cell(t_add, "series"),
      "P", fmt_p(t_add$p[t_add$term == "series"]), "\n")
  cat("Global test of the 6 interactions: CR2 Wald P", fmt_p(p_global_cr2),
      "| likelihood-ratio P", fmt_p(p_global_lr), "\n")

  rows <- data.frame(
    Term = c(sprintf("%s: case reports n = %d (%d reports), case series n = %d (%d reports)",
                     title, n_cr, k_cr, n_cs, k_cs),
             unname(SHORT[COVARS]), "Case series vs case report"),
    Main = c("", vapply(COVARS, function(v) tr_cell(t_main, v), ""), "—"),
    Added = c("", vapply(COVARS, function(v) tr_cell(t_add, v), ""), tr_cell(t_add, "series")),
    CR = c("", vapply(COVARS, function(v) tr_cell(strat[[1]], v), ""), "—"),
    CS = c("", vapply(COVARS, function(v) tr_cell(strat[[2]], v), ""), "—"),
    P_int = c("", vapply(int_terms, function(v) fmt_p(t_int$p[t_int$term == v]), ""),
              fmt_p(t_add$p[t_add$term == "series"])),
    stringsAsFactors = FALSE)
  list(rows = rows, p_cr2 = p_global_cr2, p_lr = p_global_lr, title = title)
}

A <- run_definition("Report_Type", "A. Primary definition (1 vs 2 or more patients)")
B <- run_definition("Report_Type_alt", "B. Robustness definition (1–2 vs 3 or more patients)")

tab <- rbind(A$rows, B$rows)
hdr <- c(Term = "Covariate",
         Main = sprintf("Main model\n(n = %d)", nrow(cc)),
         Added = "Report type\nadded",
         CR = "Case reports\nonly",
         CS = "Case series\nonly",
         P_int = "P")
out <- tab; names(out) <- hdr
write.csv(out, file.path(OUT_TAB, "TableS9_report_type.csv"), row.names = FALSE)

panel_rows <- which(tab$Main == "" & tab$Added == "")
foot <- c(
  "Values are time ratios (TR) with 95% CIs from the log-normal accelerated failure time model; CIs use CR2 cluster-robust standard errors clustered by source report.",
  "The last row of each panel gives the TR for case series vs case reports when report type is added to the main model.",
  "P: for each covariate, CR2 Wald test of the covariate × report type interaction in a model with all six interactions; for case series vs case report, CR2 Wald test of report type in the main model with report type added.",
  sprintf("Global test of all six interactions (CR2 Wald): primary definition P = %s; robustness definition P = %s.",
          fmt_p(A$p_cr2), fmt_p(B$p_cr2)),
  "A dash means the term is not part of that model.")
save_docx_table(out, "Table S9. Sensitivity analysis: case reports vs case series", foot,
                file.path(OUT_TAB, "TableS9_report_type.docx"),
                header_rows = panel_rows,
                level_rows = setdiff(seq_len(nrow(tab)), panel_rows))
print(out, row.names = FALSE)
cat("\nSaved Table S9\n")
