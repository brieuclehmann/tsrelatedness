#!/bin/bash -l
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=12:00:00
#$ -l mem=8G
#$ -t 1-22
#$ -o logs/balsac_grm.out
#$ -e logs/balsac_grm.err
#$ -M b.lehmann@ucl.ac.uk
#$ -m beas

# load software and environments
module load python/3.9
source .venv/bin/activate

echo "Job started at: `date`"

python3 code/balsac/07_compute_grm.py \
-chr $SGE_TASK_ID \
-rep 1
#-censor

echo "Job finished with exit code $? at: `date`"
