#!/bin/bash -l
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=12:00:00
#$ -l mem=4G
#$ -t 1-55
#$ -o logs/pca_benchmarking.out
#$ -e logs/pca_benchmarking.err
#$ -M b.lehmann@ucl.ac.uk
#$ -m beas

# load software and environments
module load python/3.9
source .venv/bin/activate

echo "Job started at: `date`"

python3 scripts/benchmarking/pca/01_base_profile.py $SGE_TASK_ID
#-censor

echo "Job finished with exit code $? at: `date`"
