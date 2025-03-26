#!/bin/bash -l

#SBATCH --job-name sample
#SBATCH --output /scratch/lanchier/2024_simulations/logs/%x_%A_%a.out
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 00:30:00
#SBATCH --array=0-215

# this script draws random samples for each replicate of each of 216 selection scenarios
# following a sampling scheme provided in a separate file

# set up modules and directories
module add r-full

main_dir=$(pwd)
work_dir=${main_dir}/data/
sim_dir="${work_dir}/simulation_Ne10000_Ngen1500"
samp_dir="${work_dir}/samples_Ne10000_Ngen1500"

# choosing the sampling scheme. this script was run with each scheme in the sample_files
# directory
sample_label="01_Ia_ideal_100g"
sampfile_dir="${main_dir}/sample_files"
sampfile="01_Ia_ideal_100g.tsv"

# creating directories for the sampliing data
mkdir -p ${samp_dir}
cd ${samp_dir}

# chose one simulation file depending on job array ID
readarray -t simfiles < ${sim_dir}/simfiles.txt
simfile=${simfiles[$SLURM_ARRAY_TASK_ID]}

# get information about values of s and initial freq. from file path
scenario_temp=${simfile#*init*/}
scenario=${scenario_temp%*/init*}
init=${scenario%*_*}
s=${scenario#*_*}

# create the corresponding directories in the sampling data folder
mkdir -p $init
cd $init
mkdir -p $s
cd $s

# run the actual sampling
Rscript ${main_dir}/supporting_scripts/sampling.R ${simfile} \
 --sampfile $sampfile_dir/$sampfile \
 --suffix $sample_label
