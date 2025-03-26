#!/bin/bash

#SBATCH --job-name gthrsim
#SBATCH --output logs/gsim_%A_%a.out
#SBATCH --error logs/gsim_%A_%a.err
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 00:20:00
#SBATCH --array=0-215

# this script gathers all 10,000 replicate of a simulated scenario into a single file

# since we make simulations for all combinations of 12 selection coefficients and 18
# initial allele frequencies, we have 12x18=216 different scenarios, hence the array
# of 216 jobs

# set up modules and directories
module add gcc r

main_dir=$(pwd)
work_dir=${main_dir}/data/simulation_Ne10000_Ngen1500

cd $work_dir

# list all 216 folders that each contain 10,000 replicates of a different simulation
folders=($(find "$PWD" -type d -name "*init*s*"))

# select one of those depending on the job array
folder=${folders[$SLURM_ARRAY_TASK_ID]}
cd $folder

# get the title of the current selection scenario
scenario=${folder#*init*/}

# paste all the replicates into a single file
paste replicates*/*.out > ${scenario}.tsv
