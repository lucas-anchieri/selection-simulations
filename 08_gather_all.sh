#!/bin/bash -l

#SBATCH --job-name all_res
#SBATCH --output logs/all_res_%j.out
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --time 00:30:00

# set label of current sample
sample_label="01_Ia_ideal_100g"

# set up modules and directories
module purge
dcsrsoft use arolle
module add gcc r

main_dir=$(pwd)
work_dir=${main_dir}/data/estimation_Ne10000_Ngen1500

cd $work_dir

# gather paths to all res and repres files
ls ./*/*/$sample_label/repres.tsv > represfiles.txt
ls ./*/*/$sample_label/res.tsv > resfiles.txt

# head and wc -l them as sanity check
head represfiles.txt
wc -l represfiles.txt
head resfiles.txt
wc -l resfiles.txt

# get all paths to repres files
readarray -t represfiles < represfiles.txt

# get the first line with col. names and create a new file with repres for all scenarios
sed '1q;d' ${represfiles[0]} > ${sample_label}_allrepres.tsv

# append all results into this file
for file in ${represfiles[@]}
do
	tail -n +2 ${file} >> ${sample_label}_allrepres.tsv
done

# get all paths to res files
readarray -t resfiles < resfiles.txt

# get the first line with col. names and create a new file with res for all scenarios
sed '1q;d' ${resfiles[0]} > ${sample_label}_allres.tsv

# append all results into this file
for file in ${resfiles[@]}
do
	tail -n +2 ${file} >> ${sample_label}_allres.tsv
done

# remove temp paths lists
rm represfiles.txt
rm resfiles.txt
