suppressPackageStartupMessages(library(data.table))

src <- Sys.getenv("SEER_ALL_HIST_FILE")
if (!nzchar(src)) stop("Set SEER_ALL_HIST_FILE to the authorised local all-histology export.")
out <- "outputs/final/supplementary"

x <- fread(
  src,
  select = c(
    "Patient ID", "Year of diagnosis", "Primary Site - labeled",
    "Histologic Type ICD-O-3", "Behavior code ICD-O-3",
    "Diagnostic Confirmation", "Age recode with single ages and 85+",
    "Survival months"
  ),
  colClasses = list(character = c("Patient ID", "Histologic Type ICD-O-3")),
  showProgress = TRUE
)

setnames(x, c("id", "year", "site", "histology", "behavior", "confirmation",
              "age_label", "months_raw"))
x[, age := suppressWarnings(as.numeric(sub("[^0-9].*$", "", age_label)))]
x[, valid_survival := grepl("^[0-9]+$", as.character(months_raw))]

audit <- x[, .(
  records = .N,
  unique_patients = uniqueN(id),
  positive_histology = sum(confirmation == "Positive histology"),
  adult_known_age = sum(!is.na(age) & age >= 18),
  valid_survival = sum(valid_survival)
), by = histology][order(-records)]

fwrite(audit, file.path(out, "expanded_comparator_histology_audit.csv"))
fwrite(x[, .N, by = site][order(site)],
       file.path(out, "expanded_comparator_site_audit.csv"))

cat("Rows:", nrow(x), "\n")
cat("Year range:", min(x$year), max(x$year), "\n")
cat("Sites:\n")
print(x[, .N, by = site][order(site)])
cat("Most frequent histologies:\n")
print(audit[1:40])
