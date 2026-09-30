b <- readRDS("work/revision_bootstrap.rds")
m <- b$replicates
p <- b$point

idx_asc <- c(1, 4, 7, 10, 13, 16)
idx_nos <- c(2, 5, 8, 11, 14, 17)
summarise_idx <- function(idx, group) data.frame(
  outcome = rep(c("OS", "CSS"), each = 3),
  landmark = rep(c(0, 12, 24), 2),
  group = group,
  estimate = p[idx],
  lower95 = apply(m[, idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
  upper95 = apply(m[, idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
)
out <- rbind(summarise_idx(idx_asc, "ASC"), summarise_idx(idx_nos, "AC_NOS"))
write.csv(out, "outputs/final/supplementary/ASC_standardized_survival_bootstrap_intervals.csv",
          row.names = FALSE)
print(out)
