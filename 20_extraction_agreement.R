# =============================================================================
# 20_extraction_agreement.R
# Agreement between the original data extraction and an independent
# re-extraction of a random sample (25 of 170 reports, one random patient each;
# sampled with seed 20261017).
#
# Run from the project root. Inputs (data/raw/reextraction/):
#   Reextraction_form_filled.xlsx                 - re-extraction by the second reviewer
#                                                   (sheet "提取表", returned by YY)
#   Reextraction_KEY_original_values_DO_NOT_SEND.xlsx - original values (JBI_Q1-Q8 pasted
#                                                   from the JBI scoring file; may be left blank)
#   Neuro_agreement_form_25_filled.xlsx           - neurological deficit judged by the second
#                                                   reviewer under the final (strict) definition
#                                                   from each patient's examination text
#                                                   (sheet "判定表（25例）"); compared with the final
#                                                   coding in data/clean/ipd_clean.csv
# Output: outputs/tables/TableS2_extraction_agreement.docx / .csv
# Packages: readxl, irr (plus those in 00_helpers.R)
#
# Categorical variables: percentage agreement and Cohen's kappa ("not reported"
# kept as its own category, so disagreement on whether a variable was reported
# counts as disagreement). Age and log(delay): two-way ICC for absolute agreement,
# single rater. JBI total score: quadratic-weighted kappa.
# =============================================================================

source(file.path("scripts", "r", "00_helpers.R"))
suppressPackageStartupMessages({ library(readxl); library(irr) })

IN <- file.path("data", "raw", "reextraction")
re  <- as.data.frame(read_excel(file.path(IN, "Reextraction_form_filled.xlsx"), sheet = "提取表"))
key <- as.data.frame(read_excel(file.path(IN, "Reextraction_KEY_original_values_DO_NOT_SEND.xlsx")))
re  <- re[!is.na(re$Patient_ID), ]
stopifnot(nrow(re) == 25, setequal(re$Patient_ID, key$Patient_ID))
re <- re[match(key$Patient_ID, re$Patient_ID), ]

# ---- recode the re-extraction to the original coding ------------------------
yn3 <- function(x) factor(c("有" = "1", "无" = "0", "无或未提及" = "0", "未报告" = "NR",
                            "仅组水平" = "NR")[x], levels = c("0", "1", "NR"))
orig3 <- function(x) factor(ifelse(is.na(x), "NR", as.character(x)), levels = c("0", "1", "NR"))
unit_to_month <- c("天" = 1 / 30, "周" = 7 / 30, "月" = 1, "年" = 12)

re_delay <- as.numeric(re[["延误：数值"]]) * unit_to_month[re[["延误：单位"]]]
re_jbi   <- as.matrix(re[, paste0("JBI Q", 1:8)])
re_jbi_total <- rowSums(re_jbi == "Yes", na.rm = TRUE)

pct <- function(a, b) mean(as.character(a) == as.character(b), na.rm = TRUE) * 100
kap <- function(a, b, w = "unweighted") {
  ok <- !is.na(a) & !is.na(b)
  if (length(unique(c(as.character(a[ok]), as.character(b[ok])))) < 2) return(NA_real_)
  kappa2(data.frame(a[ok], b[ok]), weight = w)$value
}
icc_ci <- function(a, b) {
  r <- icc(data.frame(a, b), model = "twoway", type = "agreement", unit = "single")
  sprintf("ICC %.2f (%.2f–%.2f)", r$value, r$lbound, r$ubound)
}
fmt_k <- function(k) ifelse(is.na(k), "Not estimable (no variation)", sprintf("κ = %.2f", k))

rows <- list()
add <- function(var, n, agree, stat) rows[[length(rows) + 1]] <<- data.frame(
  Variable = var, `Pairs, n` = n, `Agreement, %` = sprintf("%.0f", agree), Statistic = stat,
  check.names = FALSE)

add("Country", 25, pct(toupper(trimws(re[["报告国家/地区（Region）"]])), toupper(key$Region)), "—")
add("Year of publication", 25, pct(as.integer(re[["发表年份（Pub_Year）"]]), key$Pub_Year), "—")
add("Age, years", 25, pct(round(as.numeric(re[["年龄（Age）"]])), round(key$Age)),
    icc_ci(as.numeric(re[["年龄（Age）"]]), key$Age))
sx <- c("男" = 1, "女" = 0)[re[["性别（Sex）"]]]
add("Sex", 25, pct(sx, key$Sex), fmt_k(kap(sx, key$Sex)))
ex <- c("有明确接触史" = 1, "无或未提及" = 0)[re[["暴露史（Exposure_Hx）"]]]
add("Animal or raw dairy exposure", 25, pct(ex, key$Exposure_Hx), fmt_k(kap(ex, key$Exposure_Hx)))
add("Diagnostic delay (log months)", sum(!is.na(re_delay)),
    # agreement: same value after conversion to months; differences of <= 0.05 months (about 1.5 days)
    # or <= 10% arise only from rounding during conversion (e.g. 10 days = 0.33 vs 0.3 months)
    mean(abs(re_delay - key$Delay_to_Dx) <= 0.05 |
         abs(re_delay - key$Delay_to_Dx) / key$Delay_to_Dx <= 0.10, na.rm = TRUE) * 100,
    icc_ci(log(re_delay), log(key$Delay_to_Dx)))
for (v in list(c("发热（Fever_Status）", "Fever_Status", "Fever during the current illness"),
               c("硬膜外脓肿（Epidural_Abscess）", "Epidural_Abscess", "Epidural abscess"),
               c("腰大肌脓肿（Psoas_Abscess）", "Psoas_Abscess", "Psoas abscess"))) {
  a <- yn3(re[[v[1]]]); b <- orig3(key[[v[2]]])
  add(v[3], 25, pct(a, b), fmt_k(kap(a, b)))
}

# Neurological deficit under the final (strict) definition: second reviewer's judgement
# vs the final coding in ipd_clean.csv (not the KEY file, which holds the earlier coding)
nf <- as.data.frame(read_excel(file.path(IN, "Neuro_agreement_form_25_filled.xlsx"),
                               sheet = "判定表（25例）"))
nf <- nf[!is.na(nf$Patient_ID), ]
stopifnot(nrow(nf) == 25, setequal(nf$Patient_ID, key$Patient_ID))
nf <- nf[match(key$Patient_ID, nf$Patient_ID), ]
ncol_j <- grep("^1\\.", names(nf))[1]                       # "1. 神经功能缺损（有/无/未报告）"
a <- yn3(trimws(nf[[ncol_j]]))
if (anyNA(a)) stop("Neuro form: answers must be 有 / 无 / 未报告")
ipd_n <- load_ipd()
b <- orig3(ipd_n$Neuro_Deficit[match(key$Patient_ID, ipd_n$Patient_ID)])
add("Neurological deficit", 25, pct(a, b), fmt_k(kap(a, b)))
# optional column: second reviewer checked the examination text against the full report
chk_j <- grep("全文核对", names(nf))
neuro_checked <- length(chk_j) == 1 && all(!is.na(nf[[chk_j]]) & nzchar(trimws(nf[[chk_j]])))
add("JBI total score (0–8)", 25, pct(re_jbi_total, key$JBI_Score),
    sprintf("Weighted κ = %.2f", kap(re_jbi_total, key$JBI_Score, "squared")))
if (all(paste0("JBI_Q", 1:8) %in% names(key)) && !all(is.na(key$JBI_Q1))) {
  for (q in 1:8) {
    a <- re_jbi[, q]; b <- key[[paste0("JBI_Q", q)]]
    add(sprintf("JBI item Q%d", q), 25, pct(a, b), fmt_k(kap(a, b)))
  }
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT_TAB, "TableS2_extraction_agreement.csv"), row.names = FALSE)

save_docx_table(
  tab,
  "Table S2. Agreement between the original data extraction and an independent re-extraction",
  c("A random sample of 25 of the 170 reports included before the update search was drawn (one randomly selected patient per report; seed 20261017), and the data were extracted again by a second reviewer blinded to the original extraction. All 25 patients remain in the final cohort.",
    "Categorical variables: percentage agreement and Cohen's kappa, with 'not reported' kept as a separate category. Age and diagnostic delay: two-way intraclass correlation coefficient for absolute agreement (single rater), delay on the log scale. Delay was re-extracted as value and unit from the source and converted to months (1 month = 30 days); values were counted as agreeing when they differed by no more than 0.05 months or 10%, a difference that arises only from rounding during conversion. JBI total score: quadratic-weighted kappa.",
    if (neuro_checked)
      "Neurological deficit (final definition: muscle weakness or cauda equina / conus medullaris syndrome on examination at diagnosis): the second reviewer checked each patient's neurological examination findings against the full report and applied the definition independently; agreement is with the final coding."
    else
      "Neurological deficit (final definition: muscle weakness or cauda equina / conus medullaris syndrome on examination at diagnosis): the second reviewer applied the definition independently to each patient's neurological examination findings, quoted verbatim with page numbers from the full report; agreement is with the final coding.",
    "Kappa cannot be estimated when one reviewer gave the same answer for all patients; percentage agreement is then the relevant measure."),
  file.path(OUT_TAB, "TableS2_extraction_agreement.docx"))
