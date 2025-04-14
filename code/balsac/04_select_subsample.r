#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(purrr)

set.seed(42)

full_df <- read_csv("data/balsac_pedigree.csv")
relative_df <- read_csv("data/balsac_relatives.csv")
founder_df <- read_csv("data/balsac_founder_depths.csv")

##########################
### Select individuals ###
##########################

region_df <- full_df |>
    select(proband = ind, proband_region)

# Sample probands so that there are at least 5 {siblings, first cousins, second cousins} from each region
relative_sample <- relative_df %>%
    left_join(region_df, by = c("proband1" = "proband")) |>
    rename(proband_region1 = proband_region) |>
    left_join(region_df, by = c("proband2" = "proband")) |>
    rename(proband_region2 = proband_region) |>
    left_join(founder_df, by = c("proband1" = "proband")) %>%
    rename(average_depth1 = average_depth) %>%
    left_join(founder_df, by = c("proband2" = "proband")) %>%
    rename(average_depth2 = average_depth) %>%
    filter(average_depth1 >= 3 & average_depth2 >= 3) %>%
    filter(proband_region1 == proband_region2) %>%
    filter(!(relationship %in% c("unclassified", "half_siblings"))) %>%
    group_by(proband_region1, relationship) %>%
    slice_sample(n = 5) %>%
    ungroup() %>%
    select(proband1, proband2) %>%
    unlist() %>%
    unique()

selected_relatives <- full_df %>%
    filter(ind %in% relative_sample) %>%
    left_join(founder_df, by = "proband") 

num_remain <- selected_relatives %>%
    group_by(proband_region) %>%
    summarise(n_sample = 50 - n())

# Continue subsampling to fill 100 individuals from each region
proband_id_list <- unique(full_df$proband)

selected_df <- full_df %>%
    left_join(founder_df, by = "proband") %>%
    filter(average_depth >= 3) %>%
#    left_join(proband_df, by = c("ind" = "proband")) %>%
    filter((generation == 0) & !is.na(proband_region)) %>%
    filter(! ind %in% relative_sample) %>%
    nest(data = -proband_region) %>%
    left_join(num_remain, by = c("proband_region")) %>%
    mutate(Sample = map2(data, n_sample, sample_n)) %>%
    unnest(Sample) %>%
    select(-c(data, n_sample)) %>%
    bind_rows(selected_relatives) %>%
    arrange(proband_region)

# Check (should be 50 distinct probands from each region)
selected_df %>%
    distinct(proband_region, ind) %>%
    group_by(proband_region) %>%
    summarise(n = n())

write_csv(selected_df, "data/balsac_subsample.csv")

# Save metadata
meta_df <- selected_df %>%
    select(proband, proband_region, average_depth)
write_csv(meta_df, "data/balsac_subsample_meta.csv")
