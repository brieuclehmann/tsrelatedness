#!/bin/bash -l
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=2:30:00
#$ -l mem=4G
#$ -t 1
#$ -o logs/balsac.out
#$ -e logs/balsac.err
#$ -M b.lehmann@ucl.ac.uk
#$ -m beas

# load software and environments
module load python/3.9
source .venv/bin/activate

dir=/home/ucakble/Projects/tsrelatedness
pedigree_name=$dir/data/balsac_pedigree.csv

echo "Job started at: `date`"

python3 code/balsac/06_simulate_trees.py \
-d $dir \
-p $pedigree_name \
-chr 3 \
-rep $SGE_TASK_ID
#-censor

echo "Job finished with exit code $? at: `date`"
