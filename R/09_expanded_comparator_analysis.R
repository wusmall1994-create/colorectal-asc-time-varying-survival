suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

src <- Sys.getenv("SEER_ALL_HIST_FILE")
if (!nzchar(src)) stop("Set SEER_ALL_HIST_FILE to the authorised local all-histology export.")
final_dir <- "outputs/final"
supp_dir <- file.path(final_dir, "supplementary")

# Colorectal adenocarcinoma morphologies used in population-based colorectal
# studies, supplemented with the newer serrated and micropapillary codes.
# Histology 8570 (adenocarcinoma with squamous metaplasia) is deliberately
# excluded to keep the comparator distinct from adenosquamous carcinoma.
adenocarcinoma_codes <- c(
  "8140", "8141", "8143", "8144", "8145", "8147",
  "8210", "8211", "8213", "8220", "8221", "8255",
  "8260", "8261", "8262", "8263", "8265",
  "8480", "8481", "8490"
)

cols <- c(
  "Patient ID", "Sequence number", "Year of diagnosis", "Sex",
  "Race recode (W, B, AI, API)", "Age recode with single ages and 85+",
  "Primary Site - labeled", "Histologic Type ICD-O-3",
  "Diagnostic Confirmation",
  "Combined Summary Stage with Expanded Regional Codes (2004+)",
  "SEER historic stage A (1973-2015)", "Survival months",
  "Vital status recode (study cutoff used)",
  "SEER cause-specific death classification",
  "Grade Recode (thru 2017)", "Grade Clinical (2018+)",
  "Grade Pathological (2018+)"
)

d <- fread(src, select = cols,
           colClasses = list(character = c("Patient ID", "Histologic Type ICD-O-3")),
           showProgress = TRUE)
setnames(d, c(
  "id", "sequence", "year", "sex", "race_raw", "age_label", "site_raw",
  "histology_code", "confirmation", "summary_stage", "historic_stage",
  "months_raw", "vital", "css_class", "grade_old", "grade_clin", "grade_path"
))

d <- d[histology_code %chin% c("8560", adenocarcinoma_codes)]
d[, age := suppressWarnings(as.numeric(sub("[^0-9].*$", "", age_label)))]
d <- d[age >= 18 & confirmation == "Positive histology" &
         grepl("^[0-9]+$", as.character(months_raw))]
d[, months := as.numeric(months_raw)]
d[months == 0, months := 0.5]
d[, sequence_order := suppressWarnings(as.integer(
  sub("^(\\d+)(st|nd|rd|th).*$", "\\1", sequence)
))]
d[sequence == "One primary only", sequence_order := 0L]
d[is.na(sequence_order), sequence_order := 99L]
d[, histology := fifelse(histology_code == "8560", "ASC", "Expanded_AC")]

# Retain the earliest eligible record per patient within histology group, then
# give ASC priority for the few patients represented in both groups, matching
# the primary analysis convention.
setorder(d, histology, id, year, sequence_order)
d <- d[, .SD[1], by = .(histology, id)]
asc_ids <- d[histology == "ASC", unique(id)]
overlap_removed <- d[histology == "Expanded_AC" & id %chin% asc_ids, .N]
d <- d[!(histology == "Expanded_AC" & id %chin% asc_ids)]

d[, os_event := as.integer(vital == "Dead")]
d[, css_event := as.integer(css_class == "Dead (attributable to this cancer dx)")]
d[, stage := {
  z <- fifelse(year < 2004, historic_stage, summary_stage)
  fifelse(grepl("^Localized", z), "Localized",
    fifelse(grepl("^Regional", z), "Regional",
      fifelse(grepl("^Distant", z), "Distant", "Unknown")))
}]
d[, grade := {
  old <- fifelse(grepl("^Well|^Moderately", grade_old), "I-II",
    fifelse(grepl("^Poorly|^Undifferentiated", grade_old), "III-IV", "Unknown"))
  recent_code <- fifelse(grade_path %chin% as.character(1:4), grade_path,
    fifelse(grade_clin %chin% as.character(1:4), grade_clin, "9"))
  recent <- fifelse(recent_code %chin% c("1", "2"), "I-II",
    fifelse(recent_code %chin% c("3", "4"), "III-IV", "Unknown"))
  fifelse(year <= 2017, old, recent)
}]
d[, subsite := {
  code <- sub("-.*$", "", site_raw)
  fifelse(code %chin% c("C18.0", "C18.2", "C18.3", "C18.4"), "Right colon",
    fifelse(code %chin% c("C18.5", "C18.6", "C18.7"), "Left colon",
      fifelse(code %chin% c("C19.9", "C20.9"), "Rectosigmoid/rectum", "Other colon")))
}]
d[, race := fifelse(race_raw == "White", "White",
  fifelse(race_raw == "Black", "Black", "Other"))]
d[, period := cut(year, breaks = c(1999, 2006, 2012, 2018),
                  labels = c("2000-2006", "2007-2012", "2013-2018"))]
d[, sequence_group := fifelse(sequence_order <= 1, "Only/first", "Later primary")]
d[, histology := relevel(factor(histology), "Expanded_AC")]
for (v in c("sex", "race", "period", "stage", "grade", "subsite", "sequence_group")) {
  set(d, j = v, value = factor(d[[v]]))
}

covars <- "ns(age, 4) + sex + race + period + stage + grade + subsite + sequence_group"
intervals <- list(c(0, 6), c(6, 12), c(12, 24), c(24, 60))
res <- list()
for (outcome in c("os_event", "css_event")) {
  for (ab in intervals) {
    a <- ab[1]; b <- ab[2]
    q <- copy(d[months > a])
    q[, interval_time := pmin(months - a, b - a)]
    q[, interval_event := as.integer(get(outcome) == 1 & months <= b)]
    fit <- coxph(as.formula(paste(
      "Surv(interval_time, interval_event) ~ histology +", covars
    )), data = q, ties = "efron")
    sm <- summary(fit)
    i <- match("histologyASC", rownames(sm$coefficients))
    res[[length(res) + 1]] <- data.table(
      outcome = outcome, interval_start = a, interval_end = b,
      expanded_ac_at_risk = sum(q$histology == "Expanded_AC"),
      asc_at_risk = sum(q$histology == "ASC"),
      expanded_ac_events = sum(q$interval_event[q$histology == "Expanded_AC"]),
      asc_events = sum(q$interval_event[q$histology == "ASC"]),
      HR = exp(sm$coefficients[i, "coef"]),
      lower95 = sm$conf.int[i, "lower .95"],
      upper95 = sm$conf.int[i, "upper .95"],
      p = sm$coefficients[i, "Pr(>|z|)"]
    )
  }
}
res <- rbindlist(res)

code_counts <- d[histology == "Expanded_AC", .N, by = histology_code][order(histology_code)]
code_counts[, proportion := N / sum(N)]
cohort_counts <- d[, .(
  patients = .N, deaths = sum(os_event), cancer_deaths = sum(css_event)
), by = histology]

fwrite(res, file.path(final_dir, "expanded_adenocarcinoma_dynamic_HR.csv"))
fwrite(res, file.path(supp_dir, "TableS13_expanded_adenocarcinoma_dynamic_HR.csv"))
fwrite(code_counts, file.path(supp_dir, "TableS13b_expanded_adenocarcinoma_codes.csv"))
fwrite(cohort_counts, file.path(supp_dir, "TableS13c_expanded_adenocarcinoma_cohort.csv"))

cat("Cross-histology comparator records removed:", overlap_removed, "\n")
print(cohort_counts)
print(code_counts)
print(res)
