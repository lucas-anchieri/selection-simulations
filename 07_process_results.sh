#!/bin/bash -l

#SBATCH --job-name consolidate
#SBATCH --output logs/consolidate_%A_%a.out
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 32G
#SBATCH --time 00:30:00
#SBATCH --array 0-59

# this launches an R script that computes summary statistics for each method in each
# scenario, as well as mode and credible interval for ABC methods, and stores everything
# in a single file per scenario

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
module add r-light
ulimit -s 16384

main_dir=$(pwd)
work_dir=${main_dir}/data/estimation_Ne10000_Ngen1500_awf/init$initval/s$sval/$sample_label/

cd $work_dir

echo $initval
echo $sval
echo $sample_label

# for reach method, run the script with the file created at the last step as input
for method in approxwf bmws slattice sr
do
	Rscript $main_dir/supporting_scripts/process_results.R \
		--resfile ${method}Est.tsv \
		--init $initval \
		--s_true $sval \
		--method $method
done
