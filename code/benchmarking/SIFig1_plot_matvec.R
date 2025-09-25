library(ggplot2)
library(readr)
library(scales)

# Read the CSV
df <- read_csv("data/matvec_timing.csv")

# Plot
p <- ggplot(df, aes(x = n / 2, y = matvec)) +
  geom_point() +
  geom_line() +
  scale_x_log10(labels = label_log()) +
  scale_y_log10(labels = label_log()) +
  labs(
    x = "Sample size (diploid)",
    y = "Execution time (seconds)"
  ) +
  geom_abline(intercept = log10(0.003), slope = 0.63, 
              linetype = "dotted", colour = "grey50", alpha = 0.5) +
  theme_minimal()

p

# Save plot to PDF
ggsave("plots/Fig4-matvec-benchmark.png", plot = p, width = 5, height = 5)
