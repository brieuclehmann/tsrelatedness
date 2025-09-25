# Repository for the manuscript 'On ARGS, pedigrees, and genetic relatedness matrices'

This repository contains all the scripts used to generate the results in the manuscript 'On ARGS, pedigrees, and genetic relatedness matrices'. There are two main sets of scripts, one used to generate the benchmarking results (see `code/benchmarking`) and one used to generate the results from the Balsac analysis (see `code/balsac`). 

## Benchmarking simulations

There are three sets of benchmarking simulations, one to compare genetic relatedness matrix (GRM) computation, one to compare principal components analysis (PCA), and one to measure the time taken for matrix-vector product calculations. For first two first run `00_set_params.py` in the corresponding folder to set up the benchmarking simulation parameters. The main script is given in `01_base_profile.py`. The bash script `base_profile.sh` contains SGE instructions to run the simulations on the UCL cluster. You will need to adjust this script if you wish to run the benchmarking simulations on your local machine or your institution's cluster.

## BALSAC analysis

The BALSAC dataset is not publicly available. However, researchers can apply for access: see (https://balsac.uqac.ca/) for information. The data analysed in this manuscript is a subset of the full BALSAC pedigree. Specifically, we used a subset of the pedigree based on the 500 individuals from five parishes (see manuscript for full details).

### Preprocessing BALSAC data

The first two scripts are for preprocessing the data provided by BALSAC. Specifically, `01_preprocess_balsac.r` renames some of the columns from the raw BALSAC data and finds the maximum genealogical depth for each individual in the pedigree, required at simulation stage. Next, `02_get_founder_depths.r` computes the average founder depth for each of the probands in the pedigree. This is used to filter out probands with an average founder depth less than three

### Computing pedigree relatedness

The next set of scripts are for computing pedigree relatedness. The script `03_compute_lange_kinship.r` computes genealogical (Lange) kinship among probands up to $k$ generations using a path-based breadth-first search. Next, `04_select_subsample.r` provides a random selection of 50 closely related and less closely related individuals from each of the five parishes (see manuscript for details). Then, the `05a_compute_prm_subsample.r` and `05b_compute_prm_full.r` scripts calculate the pedigree relatedness matrix for the subsample of 250 individuals and for the full set of probands, respectively. Note that the latter script is reasonably computationally intensive.

### Simulating tree sequences

The script `06_simulate_trees.py` simulates tree sequences from a fixed pedigree. This is lightly adapted from the [genome simulator](https://github.com/LukeAndersonTrocme/genome_simulations) toolbox.  The script takes as input a text pedigree then simulates (and saves) two tree sequences: a non-recapitated one and a recapitated one. This can be repeated using different random seeds (via the --rep argument) to generate multiple tree sequences from the same pedigree. The corresponding `sge06_simulate_trees.sh` and `sge06_simulate_trees_by_chr.sh` are SGE cluster scripts to simulate 100 tree sequences for chromosome 3, with different random seeds, and a single tree sequence for each chromosome, respectively. 

### Computing genetic relatedness

Given the tree sequences, `07_compute_grm.py` computes the branch-based genetic relatedness matrix (GRM) using `tskit`. The corresponding `sge07a_compute_grm.sh` and `sge07b_compute_grm_by_chr.sh` are SGE cluster scripts to compute GRM across the 100 simulated tree sequences for chromosome 3 and across the whole genome respectively. The helper script `07a_combine_sims.py` combines the output from the GRMs across the chromosome 3 simulations, while `07b_combine_grm.py` computes the whole-genome GRM based on the per-chromosome GRMs. 

### Performing principal component analysis (PCA)

The last pair of analysis scripts performs branch PCA directly on the simulated tree sequeunces with all of the probands. The script `08_run_pca.py` runs PCA for chromosome 3 using `ts.pca`, while `08b_run_pca_genome.py` runs branch PCA on the whole genome using `ts.genetic_relatedness_vector` as a linear operator to perform eigendecomposition within `scipy.sparse.linalg.eigsh`.

### Plotting scripts

The scripts `Fig*.py` are included to reproduce Figures 1-3 in the main manuscripts.

## Citations

The genome simulations are based on the work described in the following papers, please cite them:

[Jerome Kelleher, Alison M Etheridge and Gilean McVean (2016), Efficient Coalescent Simulation and Genealogical Analysis for Large Sample Sizes, PLOS Comput Biol 12(5): e1004842. doi: 10.1371/journal.pcbi.1004842](http://dx.doi.org/10.1371/journal.pcbi.1004842)

[Dominic Nelson, Jerome Kelleher, Aaron P. Ragsdale, Claudia Moreau, Gil McVean and Simon Gravel (2020), Accounting for long-range correlations in genome-wide simulations of large cohorts, PLOS Genetics 16(5): e1008619. https://doi.org/10.1371/journal.pgen.1008619](https://doi.org/10.1371/journal.pgen.1008619)

[Luke Anderson-Trocmé,  Dominic Nelson, Shadi Zabad, Alex Diaz-Papkovich, Ivan Kryukov, Nikolas Baya, Mathilde Touvier, Ben Jeffery, Christian Dina, Hélène Vézina, Jerome Kelleher, and Simon Gravel (2023), On the genes, genealogies, and geographies of Quebec. Science 380,849-855. DOI:10.1126/science.add5300](https://doi.org/10.1126/science.add5300)