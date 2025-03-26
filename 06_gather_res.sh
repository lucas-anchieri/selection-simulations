#!/bin/bash -l

#SBATCH --job-name gather_res
#SBATCH --output logs/gather_res_%A_%a.out
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 64G
#SBATCH --time 03:30:00
#SBATCH --array 0-59

# This script lists all the result files for all replicates and all methods for a
# given scenario. Here, we consider the 12 selection coefficients and only 5 of the
# initial allele frequencies, hence the need of an array of 60 jobs to cover all
# combinations

# set up the values for initial frequency and selection coefficient
initvals=(0.01 0.05 0.1 0.25 0.5)
svals=(0 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)

# assign each job in the array to a specific combination
initid=$(($SLURM_ARRAY_TASK_ID/12))
sid=$(($SLURM_ARRAY_TASK_ID%12))

initval=${initvals[$initid]}
sval=${svals[$sid]}

#initval=$1
#sval=$2
#sample_label=$3

# initval=0.1
# sval=0.02
# sample_label="01_Ia_ideal_100g"


# set up modules and directories
module purge
dcsrsoft use arolle
module add gcc r

main_dir=$(pwd)
work_dir=${main_dir}/data/estimation_Ne10000_Ngen1500_awf/init$initval/s$sval/$sample_label/

cd $work_dir

echo $initval
echo $sval
echo $sample_label

# gather resfiles for each method
find "$PWD" -name "*approxwf.out" > approxwf_files.txt
find "$PWD" -name "*slattice.out" > slattice_files.txt
find "$PWD" -name "*sr.out.param.gz" > sr_files.txt
find "$PWD" -name "*bmws.out" > bmws_files.txt

# gather all replicates of the results for this specific scenario in a single file
Rscript $main_dir/supporting_scripts/gather_results.R
