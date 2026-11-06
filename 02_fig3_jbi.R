# =============================================================================
# 02_fig3_jbi.R
# Fig. 3. JBI critical appraisal (case report checklist), item by item
#
# Run from the project root (the folder that contains data/ and scripts/).
# Input : data/raw/jbi/*.xlsx        every JBI scoring workbook (sheet "Scoring")
#         data/clean/ipd_clean.csv   analysis cohort (Patient_ID, JBI_Score, Report_Type)
# Output: data/clean/jbi_items.csv                 item-level ratings, one row per patient
#         outputs/tables/TableS3_JBI_items.csv      n (%) per item, overall and by report type
#         outputs/figures/Fig3_JBI.pdf / .png / .tiff
# Packages: readxl, ggplot2
# =============================================================================

library(readxl)
library(ggplot2)

JBI_DIR    <- file.path("data", "raw", "jbi")
CLEAN_FILE <- file.path("data", "clean", "ipd_clean.csv")

ITEMS <- c(
  Q1 = "Q1. Demographic characteristics",
  Q2 = "Q2. History and timeline",
  Q3 = "Q3. Clinical condition at presentation",
  Q4 = "Q4. Diagnostic tests and results",
  Q5 = "Q5. Intervention or treatment",
  Q6 = "Q6. Post-intervention condition",
  Q7 = "Q7. Adverse or unanticipated events",
  Q8 = "Q8. Takeaway lessons"
)
RATINGS <- c("Yes", "Unclear", "No", "Not applicable")
# Colours checked for colour-vision deficiency; percentages are also printed on the bars
FILL <- c("Yes" = "#2E7D4F", "Unclear" = "#B8860B", "No" = "#B03A2E",
          "Not applicable" = "#8C8C8C")
# The template's yellow worked example row (not a study patient)
EXAMPLE_COMMENT <- "Post-intervention course described briefly"

# -----------------------------------------------------------------------------
# 1. Read every scoring workbook
# -----------------------------------------------------------------------------
files <- list.files(JBI_DIR, pattern = "\\.xlsx$", full.names = TRUE)
files <- files[!grepl("^~\\$", basename(files))]          # skip Excel lock files
if (length(files) == 0) stop("No .xlsx files found in ", JBI_DIR)

read_scoring <- function(f) {
  d <- read_excel(f, sheet = "Scoring", col_types = "text")
  if (ncol(d) < 13) stop(basename(f), ": expected 13 columns in sheet 'Scoring'")
  d <- d[, 1:13]
  names(d) <- c("No", "Study_ID", "Patient_ID", names(ITEMS), "Overall", "Comments")
  d <- d[!is.na(d$Patient_ID) & trimws(d$Patient_ID) != "", ]
  # Drop the template example row (first data row with the template's comment)
  if (nrow(d) > 0 && !is.na(d$Comments[1]) && startsWith(d$Comments[1], EXAMPLE_COMMENT)) {
    cat(basename(f), ": template example row removed (", d$Patient_ID[1], ")\n", sep = "")
    d <- d[-1, ]
  }
  d$source_file <- basename(f)
  cat(basename(f), ":", nrow(d), "patients\n")
  d
}
jbi <- do.call(rbind, lapply(files, read_scoring))
jbi$Patient_ID <- trimws(jbi$Patient_ID)

# Ratings: tolerate case and stray spaces, refuse anything else
norm_rating <- function(x) {
  y <- tolower(trimws(x))
  out <- ifelse(y == "yes", "Yes",
         ifelse(y == "no", "No",
         ifelse(y == "unclear", "Unclear",
         ifelse(y %in% c("na", "n/a", "not applicable"), "Not applicable", NA))))
  out
}
for (q in names(ITEMS)) {
  raw <- jbi[[q]]
  jbi[[q]] <- norm_rating(raw)
  bad <- is.na(jbi[[q]])
  if (any(bad)) {
    stop(q, ": missing or unrecognised ratings for ",
         paste0(jbi$Patient_ID[bad], " [", raw[bad], "]", collapse = ", "))
  }
}

# -----------------------------------------------------------------------------
# 2. Checks against the analysis cohort
# -----------------------------------------------------------------------------
dup <- unique(jbi$Patient_ID[duplicated(jbi$Patient_ID)])
if (length(dup)) stop("Patient_ID appears more than once in the JBI files: ",
                      paste(dup, collapse = ", "))

ipd <- read.csv(CLEAN_FILE, na.strings = "NA", stringsAsFactors = FALSE)

# Patients removed at case level after JBI scoring (sheet "Excluded_cases" of the raw
# data file) are dropped from the JBI ratings; any other mismatch still stops the script
RAW_FILE <- file.path("data", "raw", "data_cleaned_final.xlsx")
if (file.exists(RAW_FILE) && "Excluded_cases" %in% excel_sheets(RAW_FILE)) {
  excl <- trimws(read_excel(RAW_FILE, sheet = "Excluded_cases", col_types = "text")$Patient_ID)
  drop <- intersect(jbi$Patient_ID, excl)
  if (length(drop)) {
    cat("Dropped from the JBI ratings (excluded at case level):", paste(drop, collapse = ", "), "\n")
    jbi <- jbi[!jbi$Patient_ID %in% drop, ]
  }
}

only_jbi <- setdiff(jbi$Patient_ID, ipd$Patient_ID)
only_ipd <- setdiff(ipd$Patient_ID, jbi$Patient_ID)
if (length(only_jbi) || length(only_ipd)) {
  stop("Patient_ID does not match ipd_clean.csv.\n",
       "  In JBI files only (", length(only_jbi), "): ", paste(only_jbi, collapse = ", "), "\n",
       "  In ipd_clean.csv only (", length(only_ipd), "): ", paste(only_ipd, collapse = ", "),
       "\nRename these in the JBI workbook so they match the analysis table exactly.")
}
cat("JBI ratings matched to all", nrow(ipd), "patients in ipd_clean.csv\n")

# JBI score recomputed from the items vs the score in the analysis table
jbi$JBI_Score_items <- rowSums(jbi[, names(ITEMS)] == "Yes")
chk <- merge(jbi[, c("Patient_ID", "JBI_Score_items")],
             ipd[, c("Patient_ID", "JBI_Score", "Report_Type")], by = "Patient_ID")
diff_score <- chk[chk$JBI_Score_items != chk$JBI_Score, ]
if (nrow(diff_score)) {
  warning(nrow(diff_score), " patients have JBI_Score in the analysis table different from ",
          "the number of Yes items:\n",
          paste0("  ", diff_score$Patient_ID, ": table ", diff_score$JBI_Score,
                 ", items ", diff_score$JBI_Score_items, collapse = "\n"),
          "\nCorrect JBI_Score in data_cleaned_final.xlsx and rerun 02_clean_data.py.")
} else {
  cat("JBI_Score in the analysis table agrees with the item ratings for every patient\n")
}
cat("JBI score from items:\n"); print(table(jbi$JBI_Score_items))

out_items <- jbi[order(jbi$Patient_ID),
                 c("Patient_ID", names(ITEMS), "JBI_Score_items", "Overall", "source_file")]
write.csv(out_items, file.path("data", "clean", "jbi_items.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Percentages per item, overall and by report type
# -----------------------------------------------------------------------------
jbi <- merge(jbi, ipd[, c("Patient_ID", "Report_Type")], by = "Patient_ID")
n_all <- nrow(jbi); n_cr <- sum(jbi$Report_Type == 1); n_cs <- sum(jbi$Report_Type == 2)
groups <- list(
  setNames(list(jbi), paste0("All patients (n = ", n_all, ")")),
  setNames(list(jbi[jbi$Report_Type == 1, ]), paste0("Case reports (n = ", n_cr, ")")),
  setNames(list(jbi[jbi$Report_Type == 2, ]), paste0("Case series (n = ", n_cs, ")"))
)
groups <- do.call(c, groups)

rows <- list()
for (g in names(groups)) {
  d <- groups[[g]]
  for (q in names(ITEMS)) {
    tab <- table(factor(d[[q]], levels = RATINGS))
    rows[[length(rows) + 1]] <- data.frame(
      group = g, item = ITEMS[[q]], rating = RATINGS,
      n = as.integer(tab), total = nrow(d), stringsAsFactors = FALSE)
  }
}
pct <- do.call(rbind, rows)
pct$pct <- 100 * pct$n / pct$total

# Supplementary table: n (%) for every item and group
tab_s <- reshape(
  transform(pct, cell = sprintf("%d (%.1f)", n, pct))[, c("group", "item", "rating", "cell")],
  idvar = c("group", "item"), timevar = "rating", direction = "wide")
names(tab_s) <- sub("^cell\\.", "", names(tab_s))
dir.create(file.path("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(tab_s, file.path("outputs", "tables", "TableS3_JBI_items.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Plot
# -----------------------------------------------------------------------------
used <- RATINGS[RATINGS %in% pct$rating[pct$n > 0]]     # drop 'Not applicable' if unused
pct  <- pct[pct$rating %in% used, ]
pct$rating <- factor(pct$rating, levels = rev(used))     # Yes drawn at the left
pct$item   <- factor(pct$item, levels = rev(ITEMS))      # Q1 at the top
pct$group  <- factor(pct$group, levels = names(groups))
pct$label  <- ifelse(pct$pct >= 8, sprintf("%.0f", pct$pct), "")

p <- ggplot(pct, aes(x = pct, y = item, fill = rating)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            colour = "white", size = 2.6, fontface = "bold") +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_manual(values = FILL[used], breaks = used, name = NULL) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(x = "Patients (%)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        strip.text         = element_text(face = "bold", size = 9),
        axis.text.y        = element_text(colour = "grey15", size = 8.5),
        legend.position    = "bottom",
        panel.spacing.x    = unit(9, "mm"),
        plot.margin        = margin(5, 12, 5, 12))

print(p)

out_dir <- file.path("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
w <- 10; h <- 4
ggsave(file.path(out_dir, "Fig3_JBI.pdf"),  p, width = w, height = h)
ggsave(file.path(out_dir, "Fig3_JBI.png"),  p, width = w, height = h, dpi = 600)
ggsave(file.path(out_dir, "Fig3_JBI.tiff"), p, width = w, height = h, dpi = 600,
       compression = "lzw")
cat("Saved to", normalizePath(out_dir), "\n")
