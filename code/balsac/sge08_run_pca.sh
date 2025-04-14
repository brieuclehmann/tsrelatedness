#!/bin/bash -l
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=12:00:00
#$ -l mem=8G
#$ -t 1
#$ -o logs/balsac_pca.out
#$ -e logs/balsac_pca.err
#$ -M b.lehmann@ucl.ac.uk
#$ -m beas

# load software and environments
module load python/3.9
source .venv/bin/activate

echo "Job started at: `date`"

python3 code/balsac/08_run_pca.py
#-censor

echo "Job finished with exit code $? at: `date`"
