suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

final_dir <- "outputs/final"
d <- fread(file.path(final_dir, "Figure1_data.csv"))
risk <- fread(file.path(final_dir, "Table2_risk_sets_events.csv"))

risk_wide <- dcast(risk, start + end ~ histology, value.var = "at_risk")
setnames(risk_wide, c("AC_NOS", "ASC"), c("n_nos", "n_asc"))
risk_wide[, interval := sprintf("%d-%d months", start, end)]
risk_wide[, interval_label := sprintf(
  "%d-%d months\nASC n=%s; NOS n=%s",
  start, end, format(n_asc, big.mark = ","), format(n_nos, big.mark = ",")
)]

d <- merge(d, risk_wide[, .(interval, interval_label)], by = "interval", all.x = TRUE)
d[, interval_label := factor(
  interval_label,
  levels = rev(risk_wide[order(start), interval_label])
)]
d[, Outcome := factor(Outcome, levels = c("Overall survival", "Cancer-specific survival"))]

palette <- c("Overall survival" = "#2A9D8F", "Cancer-specific survival" = "#D96C5F")
p <- ggplot(d, aes(HR, interval_label, colour = Outcome)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95),
                 height = 0.13, position = position_dodge(width = 0.38), linewidth = 0.65) +
  geom_point(size = 2.4, position = position_dodge(width = 0.38)) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
                limits = c(0.5, 3), labels = c("0.5", "0.75", "1", "1.5", "2", "3")) +
  scale_colour_manual(values = palette) +
  labs(
    x = "Adjusted hazard ratio (ASC vs adenocarcinoma NOS; log scale)",
    y = NULL, colour = NULL
  ) +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.y = element_text(colour = "black", lineheight = 0.95),
    axis.text.x = element_text(colour = "black"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(file.path(final_dir, "Figure1_time_varying_HR.png"), p,
       width = 7.2, height = 4.5, dpi = 600, bg = "white")
ggsave(file.path(final_dir, "Figure1_time_varying_HR.tiff"), p,
       width = 7.2, height = 4.5, dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(final_dir, "Figure1_time_varying_HR.svg"), p,
       width = 7.2, height = 4.5, bg = "white")

fwrite(d, file.path(final_dir, "Figure1_data.csv"))
