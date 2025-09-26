library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(patchwork)

#############
#### GRM ####
#############

grm_file <- "output/benchmarking/grm_results.csv"
n_sim <- 20
if (file.exists(grm_file)) {
  sim_df <- read_csv(grm_file)
} else {
  for (i in 1:n_sim) {
    this_file <- paste0("output/benchmarking/grm/sim", i-1, ".csv")
    if (!file.exists(this_file)) {
      next
    }
    temp_df <- read_csv(this_file)
    if (i == 1) {
      sim_df <- temp_df
    } else {
      sim_df <- rbind(sim_df, temp_df)
    }
  }
  write_csv(sim_df, grm_file)
}


# Default values
default_sample_size <- 2^10
default_sequence_length <- 10^7
default_mutation_rate <- 10^-8

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colours <- gg_color_hue(3)

sim_df <- sim_df %>%
  filter(method != "ts.summary_stat")

# PLOT NUMBER OF SAMPLES
sim_df$Method <- sim_df$method
plot_df <- sim_df %>%
  filter(seq_length == default_sequence_length) %>%
  group_by(Method, num_samples) %>%
  summarise(
    time_min = min(time),
    time_max = max(time),
    time = mean(time),
    .groups = 'drop'
  )

rightmost_labels_p1 <- plot_df %>%
  group_by(Method) %>%
  filter(num_samples == max(num_samples))

p1 <- ggplot(plot_df, aes(x = num_samples, y = time, colour = Method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = time_min, ymax = time_max), width = 0.05) +
  scale_x_log10(breaks = 10^seq(1, 6), labels = label_log(), 
                limits = c(100,10000)) + 
  scale_y_log10(labels = label_log()) +
  labs(x = "Number of Samples", y = "Time (s)") +
  theme_minimal() +
  # Add reference lines
  geom_abline(intercept = log10(0.00032), slope = 2, 
              linetype = "dotted", colour = "grey50", alpha = 0.5)

p1 <- p1 +
  geom_text_repel(
    data = rightmost_labels_p1,
    aes(label = scales::comma(time, accuracy = 0.1)),  # or use 
    hjust = -0.2,  # moves the text slightly to the right
    vjust = 0.5,
    size = 3,
    show.legend = FALSE
  ) +
  coord_cartesian(clip = "off") 

# PLOT SEQUENCE LENGTH
plot_df <- sim_df %>%
  filter(num_samples == default_sample_size) %>%
  group_by(Method, seq_length) %>%
  summarise(
    time_min = min(time),
    time_max = max(time),
    time = mean(time),
    .groups = 'drop'
  )

rightmost_labels_p2 <- plot_df %>%
  group_by(Method) %>%
  filter(seq_length == max(seq_length))

p2 <- ggplot(plot_df, aes(x = seq_length, y = time, colour = Method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = time_min, ymax = time_max), width = 0.08) +
  scale_x_continuous(trans = "log10",
                     breaks = 10^seq(4, 8),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = label_log()) +
  labs(x = "Sequence length", y = "Time (s)") +
  theme_minimal()  +
  # Add reference lines
  geom_abline(intercept = log10(0.00004), slope = 1, 
              linetype = "dotted", colour = "grey50", alpha = 0.5)

p2 <- p2 +
  geom_text(
    data = rightmost_labels_p2,
    aes(label = scales::comma(time, accuracy = 0.1)),
    hjust = -0.2,
    vjust = 0.5,
    size = 3,
    show.legend = FALSE
  ) +
  coord_cartesian(clip = "off")

p3 <- p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        plot.margin = margin(1, 15, 1, 1, unit = "pt")) &
  plot_annotation(tag_levels = 'A') &
  scale_color_brewer(palette = "Set2", type = "qual")

ggsave(filename = paste0("plots/SIFig_benchmarking_plot.png"), plot = p3, width = 9, height = 5)

#############
#### PCA ####
#############

pca_file <- "output/benchmarking/pca_results.csv"

n_sim <- 20
if (file.exists(pca_file)) {
  sim_df <- read_csv(pca_file)
} else {
  for (i in 1:n_sim) {
    this_file <- paste0("output/benchmarking/pca/sim", i-1, ".csv")
    if (!file.exists(this_file)) {
      next
    }
    temp_df <- read_csv(this_file)
    if (i == 1) {
      sim_df <- temp_df
    } else {
      sim_df <- rbind(sim_df, temp_df)
    }
  }
  write_csv(sim_df, pca_file)
}

# Default values
default_sample_size <- 2^10
default_sequence_length <- 10^7
default_mutation_rate <- 10^-8

# PLOT NUMBER OF SAMPLES
plot_df <- sim_df %>%
  filter(seq_length == default_sequence_length, mutation_rate == default_mutation_rate) %>%
  group_by(method,mode,  num_samples) %>%
  summarise(
    time_min = min(time),
    time_max = max(time),
    time = mean(time),
    .groups = 'drop'
  )

max_samples = max(plot_df$num_samples)
rightmost_labels_p4 <- plot_df %>%
  group_by(method, mode) %>%
  filter(num_samples == max(num_samples))

p4 <- ggplot(plot_df, aes(x = num_samples, y = time, 
                          color = method, linetype = mode)) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_errorbar(aes(ymin = time_min, ymax = time_max), width = 0.1) +
  scale_x_continuous(trans = "log10",
                     breaks = 10^seq(1, 8),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = label_log()) +
  labs(x = "Number of Samples", y = "Time (s)") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  # Add reference lines
  geom_abline(intercept = log10(0.00032), slope = 1, 
              linetype = "dotted", colour = "grey50", alpha = 0.5)

p4

p4 <- p4 +
  geom_text(
    data = rightmost_labels_p4,
    aes(label = scales::comma(time, accuracy = 0.1)),
    hjust = -0.2,
    vjust = 0.5,
    size = 3,
    show.legend = FALSE
  ) +
  coord_cartesian(clip = "off")


# PLOT SEQUENCE LENGTH
plot_df <- sim_df %>%
  filter(num_samples == default_sample_size, mutation_rate == default_mutation_rate) %>%
  group_by(method, mode, seq_length) %>%
  summarise(
    time_min = min(time),
    time_max = max(time),
    time = mean(time),
    .groups = 'drop'
  )

rightmost_labels_p5 <- plot_df %>%
  group_by(method, mode) %>%
  filter(seq_length == max(seq_length))

p5 <- ggplot(plot_df, aes(x = seq_length, y = time, 
                          color = method, linetype = mode)) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_errorbar(aes(ymin = time_min, ymax = time_max), width = 0.1) +
  scale_x_continuous(trans = "log10",
                     breaks = 10^seq(4, 8),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = label_log()) +
  labs(x = "Sequence length", y = "Time (s)") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  # Add reference lines
  geom_abline(intercept = log10(0.0000001), slope = 1, 
              linetype = "dotted", colour = "grey50", alpha = 0.5)


p5 <- p5 +
  geom_text_repel(
    data = rightmost_labels_p5,
    aes(label = scales::comma(time, accuracy = 0.1)),
    hjust = 0.2,
    nudge_y = 0.2,
    size = 3,
    show.legend = FALSE,
    min.segment.length = 5
  ) +
  coord_cartesian(clip = "off")

p6 <- p4 + p5 + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        plot.margin = margin(1, 15, 1, 1, unit = "pt")) &
  plot_annotation(tag_levels = 'A') &
  scale_color_brewer(palette = "Set1", type = "qual")


p7 <- p3 / p6 + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = margin(1, 15, 1, 1, unit = "pt"))

ggsave(filename = paste0("plots/Fig2_benchmarking_plot.png"), plot = p6, width = 9, height = 5)
