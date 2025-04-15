library(readr)
library(dplyr)
library(ggplot2)
library(scales)

set.seed(42)

theme_set(theme_bw())

full_df <- read_csv("path_to_combined_grms")

full_df <- full_df |>
  mutate(relationship = if_else(proband1 == proband2, 'self', relationship)) |>
  mutate(relationship = if_else(is.na(relationship) , "unclassified", relationship))

full_df$relationship <- factor(full_df$relationship, 
                               levels = rev(c("self", "full_siblings", "half_siblings", "first_cousins", 
                                          "second_cousins", "3_degree_relatives", "unclassified")),
                               labels = rev(c("self", "full siblings", "half siblings", "first cousins", 
                                              "second cousins", "third cousins", "other")))
sub_df <- full_df |>
  filter(relationship != "other") |>
  filter(type == "branch_recap" & (proband_region1 == proband_region2)) |>
  group_by(relationship, proband_region1) |>
  slice_sample(prop = 0.05) |>
  ungroup()

na_df <- full_df |>
  filter(relationship == "other") |>
#  filter(prm > 1e-10 & is.na(generation)) |>
  filter(type == "branch_recap" & (proband_region1 == proband_region2)) |>
  group_by(proband_region1) |>
  slice_sample(prop = 0.002) |>
  ungroup()

sub_df <- bind_rows(sub_df, na_df)

overall_diversity <- 26367.85
theory_df <- tibble(prm = 2^seq(0, -15, by = -0.1)) |>
  mutate(grm = overall_diversity * prm / 2)

p1<- sub_df |>
  ggplot(aes(x=prm,
             y=grm, 
             group=interaction(proband1,proband2), 
             colour=relationship)) +
  geom_line(aes(x=prm,y=grm), data=theory_df, inherit.aes = FALSE, 
            color = 'blue', linetype='dotted', linewidth = 0.5) +
  geom_boxplot(width=0.08, outlier.size = 0.1, linewidth=0.3,
               position = position_dodge(width=0)) +
  scale_x_continuous(trans = "log2",
                     breaks = 2^seq(0, -15, by = -1),
                     labels = trans_format("log2", math_format(2^.x)),
                     limits = c(2^-11, 1.15)) +
  theme(legend.position="bottom", 
        panel.grid.minor = element_blank()) +
  labs(x = "Pedigree relatedness", y = "Branch relatedness") +
  guides(colour = guide_legend(nrow = 1, title = "Relationship")) +
  scale_colour_viridis_d()

             
ggsave("plots/branch_recap_sim_boxplot_combined_behind.pdf", p1, height = 4, width = 7, dpi=500)

