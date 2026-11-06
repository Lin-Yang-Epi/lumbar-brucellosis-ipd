# =============================================================================
# 00_helpers.R
# Shared settings and functions for the sensitivity analyses (07_ onwards).
# Each script sources this file; run scripts from the project root.
# Packages: clubSandwich, flextable
# =============================================================================

suppressPackageStartupMessages({
  library(clubSandwich)
  library(flextable)
})

N_EXPECTED <- 249

# Pre-specified covariates of the main model and their short labels
COVARS <- c("age10", "Sex", "Exposure_Hx", "Endemic_current", "Imaging_Era", "Fever_Status")
SHORT <- c(age10 = "Age (per 10 years)", Sex = "Male vs female",
           Exposure_Hx = "Exposure vs no exposure",
           Endemic_current = "Endemic vs non-endemic",
           Imaging_Era = "Published 2000 or later vs before",
           Fever_Status = "Fever vs no fever")

OUT_TAB <- file.path("outputs", "tables")
OUT_FIG <- file.path("outputs", "figures")
OUT_MOD <- file.path("outputs", "models")
for (d in c(OUT_TAB, OUT_FIG, OUT_MOD)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

load_ipd <- function() {
  ipd <- read.csv(file.path("data", "clean", "ipd_clean.csv"),
                  na.strings = "NA", stringsAsFactors = FALSE)
  if (nrow(ipd) != N_EXPECTED) stop("ipd_clean.csv has ", nrow(ipd), " rows; expected ", N_EXPECTED)
  ipd$age10 <- ipd$Age / 10
  ipd
}

# CR2 cluster-robust results (Satterthwaite df), one row per coefficient
robust_tab <- function(fit, cluster) {
  ct <- as.data.frame(coef_test(fit, vcov = "CR2", cluster = cluster, test = "Satterthwaite"))
  df <- ct[[grep("^df", names(ct))[1]]]
  p  <- ct[[grep("^p_", names(ct))[1]]]
  data.frame(term = names(coef(fit)), beta = ct$beta, se = ct$SE, df = df, p = p,
             lo = ct$beta - qt(0.975, df) * ct$SE,
             hi = ct$beta + qt(0.975, df) * ct$SE, stringsAsFactors = FALSE)
}

# CR2 Wald test (HTZ small-sample correction) that a set of coefficients is zero
wald_p <- function(fit, terms, cluster) {
  wt <- as.data.frame(Wald_test(fit, constraints = constrain_zero(terms),
                                vcov = "CR2", cluster = cluster, test = "HTZ"))
  wt[[grep("^p", names(wt))[1]]]
}

# Harrell's restricted cubic spline, non-linear part only (the linear term is added separately)
rcs_basis <- function(x, knots, prefix = "nl") {
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

fmt_tr <- function(b, lo, hi) sprintf("%.2f (%.2f–%.2f)", exp(b), exp(lo), exp(hi))
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# TR (95% CI) for one term of a robust_tab() result; dash if the term is absent
tr_cell <- function(tab, term) {
  r <- tab[tab$term == term, ]
  if (nrow(r) == 0) "—" else fmt_tr(r$beta, r$lo, r$hi)
}

# Covariates that vary in a data subset (a constant covariate cannot be estimated)
varying <- function(d, vars) vars[vapply(vars, function(v) length(unique(d[[v]])) > 1, logical(1))]

# Main-model formula, optionally with extra terms
main_formula <- function(covars = COVARS, extra = NULL) reformulate(c(covars, extra), "Log_Delay")

# Word table on a landscape A4 page, fixed to the page width.
#   header_rows: panel headings, bold and merged across all columns
#   level_rows : indented rows
#   bold_rows  : rows whose first-column label is bold (top-level rows)
#   first_col  : width of the first column (inches); P-value columns are narrow,
#                the remaining width is shared equally by the other columns.
# Use "\n" in column names to break long headers over two lines.
PAGE_WIDTH <- 10.2   # A4 landscape (11.69 in) minus 0.75 in margins
save_docx_table <- function(df, caption, foot, file, header_rows = integer(0),
                            level_rows = integer(0), first_col = 2.3, widths = NULL,
                            bold_rows = integer(0)) {
  nc <- ncol(df)
  labels <- names(df)

  # Coding note for the history and symptom variables: added automatically to any
  # table that shows "no exposure" or "no fever" (as a row or a column label), unless
  # the footnote already explains it
  txt <- c(labels, as.character(df[[1]]))
  hx <- c(exposure = any(grepl("no exposure", txt, ignore.case = TRUE)),
          fever    = any(grepl("no fever",    txt, ignore.case = TRUE)))
  if (any(hx) && !any(grepl("not mentioned", foot, ignore.case = TRUE))) {
    what <- paste(c("No exposure", "no fever")[hx], collapse = " / ")
    what <- paste0(toupper(substr(what, 1, 1)), substring(what, 2))
    foot <- c(foot, paste0(what, ": not mentioned or explicitly absent in the report."))
  }

  if (is.null(widths)) {
    is_p <- grepl("^P( |$|\n)", labels)
    w <- rep(NA_real_, nc); w[1] <- first_col; w[is_p] <- 0.6
    rest <- is.na(w); w[rest] <- (PAGE_WIDTH - sum(w, na.rm = TRUE)) / sum(rest)
  } else {
    if (length(widths) != nc) stop("widths needs one value per column (", nc, ")")
    w <- widths
  }

  # flextable needs unique column keys; headers such as "P" may repeat,
  # so use internal keys and show the original names as header labels
  keys <- paste0("col", seq_len(nc))
  names(df) <- keys
  ft <- flextable(df) |>
    set_header_labels(values = setNames(as.list(labels), keys)) |>
    set_caption(caption = caption) |>
    bold(part = "header") |>
    align(j = seq(2, nc), align = "center", part = "all") |>
    valign(valign = "bottom", part = "header") |>
    add_footer_lines(foot) |>
    font(fontname = "Times New Roman", part = "all") |>
    fontsize(size = 9, part = "all") |>
    padding(padding.top = 1, padding.bottom = 1, part = "body")
  if (length(level_rows)) ft <- padding(ft, i = level_rows, j = 1, padding.left = 12, part = "body")
  if (length(bold_rows)) ft <- bold(ft, i = bold_rows, j = 1, part = "body")   # top-level rows
  if (length(header_rows)) {
    ft <- merge_h_range(ft, i = header_rows, j1 = 1, j2 = nc, part = "body")
    ft <- bold(ft, i = header_rows, part = "body")
    ft <- align(ft, i = header_rows, align = "left", part = "body")
  }
  ft <- width(ft, j = seq_len(nc), width = w)
  ft <- set_table_properties(ft, layout = "fixed")
  saveRDS(ft, sub("\\.docx$", ".rds", file))   # read by 23_supplement.R
  sect <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(left = 0.75, right = 0.75, top = 0.75, bottom = 0.75))
  save_as_docx(ft, path = file, pr_section = sect)
  invisible(ft)
}
