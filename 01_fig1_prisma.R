# =============================================================================
# 01_fig1_prisma.R
# Fig. 1. PRISMA-IPD flow diagram
#
# Run from the project root (the folder that contains data/ and scripts/).
# Input : data/clean/ipd_clean.csv   (used only to check n = 249 and to count
#                                     the studies contributing to the cohort)
#
# The original search (to 22 May 2026) and the update search (19 Oct 2026, records added
# to the databases since the original search) are combined in one flow; each box also
# shows the update-search count so both searches can be traced.
# Output: outputs/figures/Fig1_PRISMA_IPD.pdf / .png / .tiff
# Packages: ggplot2 only
# =============================================================================

library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Counts (from the screening records; edit here only)
# -----------------------------------------------------------------------------
# Original search (to 22 May 2026)
o <- list(
  db_pubmed = 1521, db_wos = 1580, db_scopus = 2937,
  other_sources   = 439,   # Zotero library of preliminary searches for related reviews, imported together
  test_set        = 68,    # articles used to test the search strategy (part of other_sources)
  duplicates      = 2973,
  already         = 0,
  screen_excluded = 2472,
  not_retrieved   = 337,   # batch 1: 195 of 699; batch 2: 142 of 333
  ft_excluded     = 456,
  citation        = 3,
  ipd_extracted   = 519    # 515 + Kim_2020, Karadag_2023, Tuz_2025, Zhang_2018
)
# Update search (19 Oct 2026; records added to the databases since the original search)
u <- list(
  db_pubmed = 33, db_wos = 31, db_scopus = 89,
  other_sources   = 0,
  duplicates      = 52,
  already         = 3,     # records already assessed in the original search
  screen_excluded = 72,
  not_retrieved   = 2,
  ft_excluded     = 18,
  citation        = 0,     # citation searching of the 5 new reports: update here if any are found
  ipd_extracted   = 6
)
for (s in c("o", "u")) {   # derived counts for each search
  x <- get(s)
  x$db_total    <- x$db_pubmed + x$db_wos + x$db_scopus
  x$screened    <- x$db_total + x$other_sources - x$duplicates - x$already
  x$sought      <- x$screened - x$screen_excluded
  x$assessed    <- x$sought - x$not_retrieved
  x$included_db <- x$assessed - x$ft_excluded
  x$studies     <- x$included_db + x$citation
  assign(s, x)
}
stopifnot(o$screened == 3504, o$sought == 1032, o$assessed == 695, o$included_db == 239,
          u$screened == 98, u$sought == 26, u$assessed == 24, u$included_db == 6)

# Combined counts used in the diagram
n <- list(
  db_total        = o$db_total + u$db_total,
  db_pubmed = o$db_pubmed, db_wos = o$db_wos, db_scopus = o$db_scopus,
  other_sources   = o$other_sources,
  test_set        = o$test_set,
  duplicates      = o$duplicates + u$duplicates,
  already         = u$already,
  screened        = o$screened + u$screened,
  screen_excluded = o$screen_excluded + u$screen_excluded,
  sought          = o$sought + u$sought,
  not_retrieved   = o$not_retrieved + u$not_retrieved,
  assessed        = o$assessed + u$assessed,
  ft_excluded     = o$ft_excluded + u$ft_excluded,
  included_db     = o$included_db + u$included_db,
  citation        = o$citation + u$citation,
  studies         = o$studies + u$studies,
  ipd_extracted   = o$ipd_extracted + u$ipd_extracted,
  no_delay        = 234,   # 231 original; + Alzahrani_2026 (update search, in-hospital timeline only)
                           # + Tuz_2025 (no symptom duration reported)
                           # + Li_2026 Case 2 (no lumbar-related symptom onset reported)
  case_excluded   = 42,
  final           = 249
)
n$case_assessed <- n$ipd_extracted - n$no_delay

ft_reasons <- c(   # full-text exclusion reasons (original + update), must sum to ft_excluded
  "No spinal infection"                                    = 237 + 10,
  "Cervical or thoracic spine only"                        = 80 + 6,
  "Aggregate data only, no individual data"                = 56,
  "Not brucellosis"                                        = 42 + 1,
  "Duplicate report"                                       = 17,
  "Non-human or archaeological material"                   = 10 + 1,
  "Review, letter or commentary"                           = 11,
  "Full text did not match the record or not traceable"    = 3
)
case_reasons <- c( # case-level exclusion reasons, must sum to case_excluded
  "Isolated sacroiliitis"                        = 32,
  "Cervical spine only"                          = 1,
  "Lumbar involvement could not be confirmed"    = 4,   # + Ishihara_2024, Jahanshahi_2024 (2026-10-21)
  "Original not retrievable (cited elsewhere)"   = 2,   # known only from a secondary citation
  "No spinal brucellosis"                        = 3   # Li_2026: vascular infection only
)

# Arithmetic checks: the script stops if any box does not add up
stopifnot(
  n$db_total + n$other_sources - n$duplicates - n$already == n$screened,
  n$screened - n$screen_excluded == n$sought,
  n$sought - n$not_retrieved == n$assessed,
  n$assessed - n$ft_excluded == n$included_db,
  n$included_db + n$citation == n$studies,
  n$ipd_extracted - n$no_delay == n$case_assessed,
  n$case_assessed - n$case_excluded == n$final,
  sum(ft_reasons) == n$ft_excluded,
  sum(case_reasons) == n$case_excluded
)

# -----------------------------------------------------------------------------
# 2. Check against the clean data and count contributing studies
# -----------------------------------------------------------------------------
clean_file <- file.path("data", "clean", "ipd_clean.csv")
if (!file.exists(clean_file)) {
  stop("Cannot find ", clean_file, ". Working directory is ", getwd(),
       ". Open the project folder (the one containing data/ and scripts/) first.")
}
ipd <- read.csv(clean_file, na.strings = "NA", stringsAsFactors = FALSE)
if (nrow(ipd) != n$final) {
  stop("ipd_clean.csv has ", nrow(ipd), " rows, but the flow diagram says ", n$final)
}
n_studies_final <- length(unique(ipd$Cluster_ID))
cat("Final cohort:", nrow(ipd), "patients from", n_studies_final, "studies\n")

# -----------------------------------------------------------------------------
# 3. Layout (plot coordinates: x 0-100, y 0-140)
# -----------------------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
reason_lines <- function(title, total, r) {
  paste0(title, " (n = ", fmt(total), "):\n",
         paste0("  ", names(r), " (n = ", fmt(r), ")", collapse = "\n"))
}

MX <- 32; MW <- 46   # main column: centre and width
RX <- 78; RW <- 40   # right column: centre and width

boxes <- data.frame(
  id = c("db", "other", "dup", "screen", "scr_ex", "sought", "notret",
         "assessed", "ft_ex", "inc_db", "citation", "studies",
         "ipd", "nodelay", "case", "case_ex", "final"),
  x  = c(MX, RX, MX, MX, RX, MX, RX, MX, RX, MX, RX, MX, MX, RX, MX, RX, MX),
  y  = c(131, 131, 117, 104, 104, 93, 93, 82, 74, 55, 55, 43, 31, 31, 20, 16, 6),
  w  = c(MW, RW, MW, MW, RW, MW, RW, MW, RW, MW, RW, MW, MW, RW, MW, RW, MW),
  h  = c(14, 10, 9, 6, 6, 6, 6, 6, 21, 6, 6, 6, 6, 6, 6, 14, 7),
  type = c("main", "main", "excl", "main", "excl", "main", "excl", "main", "excl",
           "main", "main", "main", "main", "excl", "main", "excl", "final"),
  align = c(rep("centre", 8), "left", rep("centre", 6), "left", "centre"),
  label = c(
    paste0("Records identified from databases (n = ", fmt(n$db_total), ")\n",
           "Original search (to 22 May 2026):\n",
           "PubMed ", fmt(o$db_pubmed), " | Web of Science ", fmt(o$db_wos),
           " | Scopus ", fmt(o$db_scopus), "\n",
           "Update search (19 Oct 2026):\n",
           "PubMed ", fmt(u$db_pubmed), " | Web of Science ", fmt(u$db_wos),
           " | Scopus ", fmt(u$db_scopus)),
    paste0("Records identified from other sources (n = ", fmt(n$other_sources), "):\n",
           "preliminary search library for related reviews\n",
           "on lumbar spinal infection, including\n",
           n$test_set, " test-set articles"),
    paste0("Records removed before screening:\nduplicate records (n = ", fmt(n$duplicates), ")\n",
           "already assessed in the original search (n = ", fmt(n$already), ")"),
    paste0("Records screened\n(n = ", fmt(n$screened), "; update search ", fmt(u$screened), ")"),
    paste0("Records excluded\n(n = ", fmt(n$screen_excluded), "; update search ",
           fmt(u$screen_excluded), ")"),
    paste0("Reports sought for retrieval\n(n = ", fmt(n$sought), "; update search ", fmt(u$sought), ")"),
    paste0("Reports not retrieved\n(n = ", fmt(n$not_retrieved), "; update search ",
           fmt(u$not_retrieved), ")"),
    paste0("Reports assessed for eligibility\n(n = ", fmt(n$assessed), "; update search ",
           fmt(u$assessed), ")"),
    paste0("Reports excluded (n = ", fmt(n$ft_excluded), "; update search ",
           fmt(u$ft_excluded), "):\n",
           paste0("  ", names(ft_reasons), " (n = ", fmt(ft_reasons), ")", collapse = "\n")),
    paste0("Studies included from databases\n(n = ", fmt(n$included_db), "; update search ",
           fmt(u$included_db), ")"),
    paste0("Studies identified by backward\ncitation searching (n = ", fmt(n$citation), ")"),
    paste0("Studies included for IPD extraction\n(n = ", fmt(n$studies), ")"),
    paste0("Individual patient data extracted\n(n = ", fmt(n$ipd_extracted), " patients)"),
    paste0("Patients excluded (n = ", fmt(n$no_delay), "):\ndiagnostic delay not reported"),
    paste0("Patients assessed for case-level eligibility\n(n = ", fmt(n$case_assessed), ")"),
    reason_lines("Patients excluded", n$case_excluded, case_reasons),
    paste0("Final analytic cohort\n(n = ", fmt(n$final), " patients from ",
           n_studies_final, " studies)")
  ),
  stringsAsFactors = FALSE
)
boxes$xmin <- boxes$x - boxes$w / 2
boxes$xmax <- boxes$x + boxes$w / 2
boxes$ymin <- boxes$y - boxes$h / 2
boxes$ymax <- boxes$y + boxes$h / 2
boxes$tx   <- ifelse(boxes$align == "left", boxes$xmin + 1.5, boxes$x)
boxes$hj   <- ifelse(boxes$align == "left", 0, 0.5)
B <- function(id) boxes[boxes$id == id, ]

# Arrows: x, y (start) -> xend, yend
down  <- function(a, b) data.frame(x = B(a)$x, y = B(a)$ymin, xend = B(b)$x, yend = B(b)$ymax)
right <- function(a, b, y = B(a)$y) data.frame(x = B(a)$xmax, y = y, xend = B(b)$xmin, yend = y)
arrows <- rbind(
  data.frame(x = B("db")$x, y = B("db")$ymin, xend = B("dup")$x, yend = B("dup")$ymax),
  data.frame(x = B("other")$x - 8, y = B("other")$ymin,
             xend = B("dup")$xmax - 4, yend = B("dup")$ymax),
  down("dup", "screen"), right("screen", "scr_ex"),
  down("screen", "sought"), right("sought", "notret"),
  down("sought", "assessed"), right("assessed", "ft_ex"),
  down("assessed", "inc_db"),
  down("inc_db", "studies"),
  data.frame(x = B("citation")$x - 8, y = B("citation")$ymin,
             xend = B("studies")$xmax - 4, yend = B("studies")$ymax),
  down("studies", "ipd"), right("ipd", "nodelay"),
  down("ipd", "case"), right("case", "case_ex"),
  down("case", "final")
)

# Phase labels on the left
phases <- data.frame(
  label = c("Identification", "Screening", "Included", "Individual patient data"),
  ymin  = c(112, 60, 39, 1.5),
  ymax  = c(138, 110, 59, 35)
)
phases$y <- (phases$ymin + phases$ymax) / 2

fill_col <- c(main = "#EAF2FB", excl = "#FBECEC", final = "#FFE699")
line_col <- "grey25"

# -----------------------------------------------------------------------------
# 4. Plot
# -----------------------------------------------------------------------------
p <- ggplot() +
  geom_rect(data = phases, aes(xmin = 1, xmax = 6, ymin = ymin, ymax = ymax),
            fill = "#D9D9D9", colour = NA) +
  geom_text(data = phases, aes(x = 3.5, y = y, label = label),
            angle = 90, size = 3, fontface = "bold") +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
               linewidth = 0.4, colour = line_col) +
  geom_rect(data = boxes[boxes$type != "final", ],
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
            colour = line_col, linewidth = 0.4) +
  geom_rect(data = boxes[boxes$type == "final", ],
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
            colour = line_col, linewidth = 0.9) +
  geom_text(data = boxes, aes(x = tx, y = y, label = label, hjust = hj),
            size = 2.75, lineheight = 0.95) +
  scale_fill_manual(values = fill_col, guide = "none") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 140), expand = FALSE) +
  theme_void()

print(p)   # shows the figure in the Plots pane

# -----------------------------------------------------------------------------
# 5. Save
# -----------------------------------------------------------------------------
out_dir <- file.path("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
w <- 7.5; h <- 10   # inches
ggsave(file.path(out_dir, "Fig1_PRISMA_IPD.pdf"),  p, width = w, height = h)
ggsave(file.path(out_dir, "Fig1_PRISMA_IPD.png"),  p, width = w, height = h, dpi = 600)
ggsave(file.path(out_dir, "Fig1_PRISMA_IPD.tiff"), p, width = w, height = h, dpi = 600,
       compression = "lzw")
cat("Saved to", normalizePath(out_dir), "\n")
