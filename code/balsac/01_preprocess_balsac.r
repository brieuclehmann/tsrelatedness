#!/usr/bin/env Rscript

library(readr)
library(dplyr)
source("utils/pedigree_tools.R")

set.seed(42)

################################
## Prepare pedigree dataframe ##
################################

pedigree_df <- read_delim("data/ascendance.csv", delim = ';')
colnames(pedigree_df) <- c('proband_region', 'proband', 'ind', 'sex', 
                        'mother', 'father', 'decade', 'ParentsMarriageLocation', 
                        'ParentsMarriageLocationID', 'Origin', 'OriginID', 'Sosa', 'Generation')

pedigree_df$mother <- as.integer(pedigree_df$mother)
pedigree_df$father <- as.integer(pedigree_df$father)

pedigree_df <- distinct(pedigree_df) # Remove small number of duplicate rows
proband_id <- unique(pedigree_df$proband)
max_generation <- maximum_genealogical_depth(pedigree_df, proband_id)

full_df <- pedigree_df %>%
    left_join(max_generation, by = "ind")
write_csv(full_df, "data/balsac_pedigree.csv")
