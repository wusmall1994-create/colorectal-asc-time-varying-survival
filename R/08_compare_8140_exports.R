suppressPackageStartupMessages(library(data.table))

cols <- c("Patient ID", "Sequence number", "Year of diagnosis",
          "Primary Site - labeled", "Histologic Type ICD-O-3",
          "Behavior code ICD-O-3", "Diagnostic Confirmation",
          "Age recode with single ages and 85+", "Survival months")

all_hist_file <- Sys.getenv("SEER_ALL_HIST_FILE")
ac_nos_file <- Sys.getenv("SEER_AC_NOS_FILE")
if (!nzchar(all_hist_file) || !nzchar(ac_nos_file)) stop("Set SEER_ALL_HIST_FILE and SEER_AC_NOS_FILE.")
new <- fread(all_hist_file,
             select = cols, colClasses = list(character = c("Patient ID", "Histologic Type ICD-O-3")),
             showProgress = FALSE)[`Histologic Type ICD-O-3` == "8140"]
old <- fread(ac_nos_file,
             select = cols, colClasses = list(character = c("Patient ID", "Histologic Type ICD-O-3")),
             showProgress = FALSE)
old <- old[`Year of diagnosis` <= 2018 &
             sub("-.*$", "", `Primary Site - labeled`) %chin%
             c("C18.0", paste0("C18.", 2:9), "C19.9", "C20.9")]

key <- c("Patient ID", "Sequence number", "Year of diagnosis", "Primary Site - labeled")
setkeyv(new, key); setkeyv(old, key)
new_only <- fsetdiff(new[, ..key], old[, ..key])
old_only <- fsetdiff(old[, ..key], new[, ..key])

report <- data.table(
  metric = c("new_8140_records", "old_8140_records_to_2018", "new_only", "old_only"),
  N = c(nrow(new), nrow(old), nrow(new_only), nrow(old_only))
)
fwrite(report, "outputs/final/supplementary/expanded_comparator_8140_export_comparison.csv")
print(report)
cat("Old-only counts by year:\n")
print(old[old_only, on = key, .N, by = `Year of diagnosis`][order(`Year of diagnosis`)])
cat("Old-only counts by site:\n")
print(old[old_only, on = key, .N, by = `Primary Site - labeled`][order(-N)])
