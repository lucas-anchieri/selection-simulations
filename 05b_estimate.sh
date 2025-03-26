#!/bin/bash -l

#SBATCH --job-name estimate
#SBATCH --output /scratch/lanchier/2024_simulations/logs/%x_%A_%a.out
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 2G
#SBATCH --time 02:30:00
#SBATCH --array=0-999

# this is the big one. For each replicate of a given scenario, runs estimation of
# selection coefficients with each of the 4 methods in turn

initval=$1
sval=$2
sample_label=$3

#initval=0.1
#sval=0.02
#sample_label="01_Ia_ideal_100g"

# setting modules and directories
module purge
dcsrsoft use arolle
module add gcc gsl slim python r

main_dir=$(pwd)
work_dir=${main_dir}/data/
samp_dir=${work_dir}/samples_Ne10000_Ngen1500
est_dir=${work_dir}/estimation_Ne10000_Ngen1500_awf
samplefile=${samp_dir}/init${initval}/s${sval}/init${initval}_s${sval}_${sample_label}.tsv

# depending on the job array id, get a different replicate every time (3rd column is the
# first replicate)
rep_fact=0
rep=$(($SLURM_ARRAY_TASK_ID + $rep_fact))
repcol=$(($rep + 3))

# create the directory the will store the results of the estimations
mkdir -p ${est_dir}
cd ${est_dir}
pwd

# get info from the path of the sample file
scenario_temp=${samplefile#*init*/s*/}
scenario=${scenario_temp%*.tsv}
IFS='_' read -r -a params <<< "$scenario"
init=${params[0]}
s=${params[1]}
printf -v sample_temp '%s_' "${params[@]:2}"
sample=${sample_temp%_}

# use this info to create necessary nested directories
mkdir -p $init
cd $init
pwd

mkdir -p $s
cd $s
pwd

mkdir -p $sample
cd $sample
pwd

# create 10 replicate folders that will each store 100 nested replicate folders
# again, we choose to nest them in order to avoid too many files being written on at the
# same time in the same directory
for i in {0..9}
do
	mkdir -p replicates$(($i * 100))-$(($i * 100 + 99))
done

# select the appropriate replicate folder for the current scenario
repfolder=replicates$((($SLURM_ARRAY_TASK_ID / 100) * 100))-$((($SLURM_ARRAY_TASK_ID/ 100) * 100 + 99))

cd $repfolder
pwd

# create the aforementioned nested replicate directory
mkdir -p ${scenario}_${SLURM_ARRAY_TASK_ID}
cd ${scenario}_${SLURM_ARRAY_TASK_ID}
pwd

# get only the current replicate from the file with the sampling data
cut -f 1,2,${repcol} ${samplefile} > sample.tsv

# format samples so they can be read by the methods
# BMWS needs a .vcf, so we start with an empty one
cp ${main_dir}/sample_files/sim.vcf ./sample_bmws.vcf
# run a script that will format everything nicely
Rscript ${main_dir}/supporting_scripts/format_samples.R sample.tsv

# apply in turn each of the method using the formatted sample files
echo approxwf
ApproxWF \
	${main_dir}/launch_methods/estimate.input \
	N=10000 \
	loci=sample_approxwf.loci \
	outName=sample_approxwf.out > sample_approxwf.log

echo slattice
slattice \
	sample_slattice.txt \
	20000 > sample_slattice.out

echo sr
sr \
	-D sample_sr.txt -G 25 -N 10000 -h 0.5 \
	-n 1200000 -d 0.001 -F 20 -f 1000 -s 100 \
	-P ${main_dir}/launch_methods/constant.pop -a \
	-o sample_sr.out > sample_sr.log

echo bmws
bmws analyze \
	--vcf sample_bmws.vcf \
	--meta sample_bmws.meta \
	-g 25 \
	-l 4.5 > sample_bmws.out

echo $repcol
