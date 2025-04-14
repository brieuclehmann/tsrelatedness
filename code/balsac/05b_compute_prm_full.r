#!/usr/bin/env Rscript

library(dplyr)
library(readr)
# library(package = "devtools")
# devtools::install_github(repo = "https://github.com/Rpedigree/pedigreeTools")
library(package = "pedigreeTools")

# LOAD DATA
full_df <- read_csv("data/balsac_pedigree.csv")

proband_df <- full_df %>%
  filter(!(ind %in% mother) & !(ind %in% father))

# Remap IDs
ped <- full_df %>%
    select(-proband) %>%
    arrange(-generation) %>%
    distinct(ind, father, mother) %>%
    mutate(id = row_number())

id_map <- ped$id
names(id_map) <- ped$ind

ped <- ped %>%
  mutate(fid = id_map[as.character(father)],
         mid = id_map[as.character(mother)]) %>%
  mutate(sample = if_else(ind %in% proband_df$ind, TRUE, FALSE))

id_map <- as.character(ped$proband)
names(id_map) <- ped$id

# ---- Calculate pedigree-base relatedness (PRM) -------------------------------

# You can skip this section and rather import exported data at the end of this section

ped2 <- pedigree(sire = ped$fid,
                 dam =  ped$mid,
                 label = ped$id)
(nInd <- length(ped2@label)) # 34_882

# Size of PRM among all individuals
nInd * nInd * 8 / 2^30 # ~9GB

# Size of PRM among samples
samplesInt <- ped$id[ped$sample]
samplesChar <- as.character(samplesInt)
(nSamp <- length(samplesInt)) # 500
nSamp * nSamp * 8 / 2^30 # ~0GB

# Build PRM among the samples
# BEWARE: you must source the version of getASubset() from
#         https://github.com/Rpedigree/pedigreeTools/pull/5
#         since sub-pedigree is not encoded 1:n
PRM_full_mat <- matrix(0, nSamp, nSamp, dimnames = list(samplesChar, samplesChar))

# Split the compute to get around memory issues
breaks <- split(samplesChar, ceiling(seq_along(samplesChar)/500))
for (i in 2:length(breaks)) {
  for (j in 1:(i-1)) {
    print(c(j,i))
    subChar <- c(unlist(breaks[as.character(i)]), unlist(breaks[as.character(j)]))
    PRM_samples <- getASubset(ped = ped2, labs = subChar)
    PRM_full_mat[subChar, subChar] <- as.matrix(PRM_samples)
  }
}

str(PRM_full_mat)
PRM_full_mat[samplesChar[1:5], samplesChar[1:5]]

# Change back to BALSAC IDs
new_names = id_map[colnames(PRM_full_mat)]
colnames(PRM_full_mat) <- new_names

write_csv(as_tibble(PRM_full_mat), "output/chr3/prm_full.csv")
