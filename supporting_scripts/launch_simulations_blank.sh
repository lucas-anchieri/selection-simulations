#!/bin/bash

#SBATCH --job-name simul
#SBATCH --output logs/simul_%A_%a.out
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 00:15:00
#SBATCH --array=0-599

# We want to run 10,000 replicates of each simulation scenario. A given batch will cover
# 6 selection coefficients. Due to limitations of the queuing system, we cannot just
# launch 60,000 jobs, so we run 100 jobs for each value of s, each containing a for loop
# running the simulation 100 times in order to recover 100x100=10,000 replicates each.
# Thus, the present script is a job array of 600 jobs, allocating 100 to each of the 6
# values of s

# load modules and get paths to the slim scripts of the current batch
module add gcc slim python r gsl
readarray -t slimfiles < slimfiles_$1.txt

# change the selection value after 100 jobs
selecID=$(($SLURM_ARRAY_TASK_ID / 100))

# replicate ID between 1 and 100 for the current selection value
repID=$(($SLURM_ARRAY_TASK_ID % 100))

# from the slim script name, get the folder where simulations are to be stored
selecfolder=${slimfiles[${selecID}]%/slim*}

# get the exact path to the slim script
simfile=${slimfiles[$selecID]#"${selecfolder}/"}

# define which replicate folder to store the data into
repfolder=replicates$((($repID / 10) * 10))-$((($repID / 10) * 10 + 9))

# create 10 temporary folders, each to store 100 temporary replicate files
scenario=${selecfolder#*init*/}
cd ${selecfolder}
mkdir ${repfolder}/${scenario}_${repID}

# run the slim simulation 100 times and store each replicate in a separate file
for i in {1..100}
do
	OUTFILE=${repfolder}/${scenario}_${repID}/slim_${scenario}_${repID}_${i}.out
	slim ${simfile} | tail -n +14 > $OUTFILE
done

# paste all 100 output files in a single file that will gather all 100 replicates
paste ${repfolder}/${scenario}_${repID}/* > ${repfolder}/slim_${scenario}_${repID}.out

# since the 100 replicates are stored in a single file, we can remove the temporary output
rm -r ${repfolder}/${scenario}_${repID}
