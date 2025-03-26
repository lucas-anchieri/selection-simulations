#!/bin/bash

#SBATCH --job-name runsim
#SBATCH --output logs/runsim_%j.out
#SBATCH --error logs/runsim_%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user lucas.anchieri@unil.ch

#SBATCH --partition cpu

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 00:30:00

# set up modules and directories
module add gcc python

main_dir=$(pwd)
work_dir=${main_dir}/data

# run the simulations for 1,500 generations and with a population size of 10,000
NGEN=1500
POP_SIZE_Ne=10000
POP_SIZE_2Ne=$(( ${POP_SIZE_Ne} * 2 ))

# We simulate allele trajectories for a range of values for s (from 0 to 0.1) and a
# range of initial allele frequencies (from 1% to 50%). Since the queueing system of
# our HPC only allows running 10,000 jobs at a time and we are simulating a lot of
# replicates, we performed the simulations in 4 batches, by running this script 4
# times with different values for s and initial frequency

## FIRST BATCH ##
BATCH=1
INIT_FREQS=("0.01" "0.02" "0.03" "0.04" "0.05" "0.06" "0.07" "0.08" "0.09")
S_VALS=("0" "0.005" "0.01" "0.02" "0.03" "0.04")

## SECOND BATCH ##
#BATCH=2
#INIT_FREQS=("0.01" "0.02" "0.03" "0.04" "0.05" "0.06" "0.07" "0.08" "0.09")
#S_VALS=("0.05" "0.06" "0.07" "0.08" "0.09" "0.1")

## THIRD BATCH ##
#BATCH=3
#INIT_FREQS=("0.1" "0.15" "0.2" "0.25" "0.3" "0.35" "0.4" "0.45" "0.5")
#S_VALS=("0" "0.005" "0.01" "0.02" "0.03" "0.04")

## FOURTH BATCH ##
# BATCH=4
# INIT_FREQS=("0.1" "0.15" "0.2" "0.25" "0.3" "0.35" "0.4" "0.45" "0.5")
# S_VALS=("0.05" "0.06" "0.07" "0.08" "0.09" "0.1")

# creating a directory where the simulations will be stored
sim_dir=${work_dir}/simulation_Ne${POP_SIZE_Ne}_Ngen${NGEN}
mkdir -p ${sim_dir} && cd $_

# for each of the initial frequency values
for INIT_FREQ in ${INIT_FREQS[@]}
do
  # compute the initial derived allele count
  INIT_tmp="$(python -c 'print('${POP_SIZE_2Ne}'*'${INIT_FREQ}')')"
  INIT=${INIT_tmp%.*}

  # create a folder for this specific initial frequency
  echo $INIT_FREQ
  scenario=init${INIT_FREQ}
  mkdir -p ${scenario}
  cd ${scenario}

  # for each value of s
  for S_VAL in ${S_VALS[@]}
  do
    # create a folder for this value of s
    echo $S_VAL
    scenario_s=${scenario}_s${S_VAL}
    mkdir -p ${scenario_s}
    cd ${scenario_s}

    # from the default slim script, create a new one and set the values for this scenario
    SLIMFILE="slim_script${BATCH}_${scenario_s}.txt"
    sed "s/S_COEFF/${S_VAL}/g" ${main_dir}/supporting_scripts/slim_default.txt > $SLIMFILE
    sed -i "s/POP_SIZE/$POP_SIZE_Ne/g" $SLIMFILE
    sed -i "s/INIT/$INIT/g" $SLIMFILE
    sed -i "s/N_GEN/$NGEN/g" $SLIMFILE

    # create 10 replicate folders
    for i in {0..9}
    do
      mkdir -p replicates$(($i * 10))-$(($i * 10 + 9))
    done
    # later on, each of these 10 folders will contain nested folder each containing part
    # of the output files of the 10,000 replicates of the simulation
    # we do this in order to avoid storing too many files in the same directory, which
    # could slow down writing of the output data
    # in the end, each of the 10 replicate folders will contain 10 files, each storing
    # data for 100 replicates, for a total of 10x10x100=10,000 replicates

    cd ..
  done

  # gather all slim scripts for all the different scenarios
  find "$PWD" -name "*slim_script${BATCH}_*.txt" > slimfiles_$BATCH.txt

  # set up the shell script that will run replicates of the slim simulations
  cp ${main_dir}/supporting_scripts/launch_simulations_blank.sh ./simulate_init${INIT_FREQ}.sh
  mkdir -p logs

  # run the simulations for the current batch
  sbatch ./simulate_init${INIT_FREQ}.sh $BATCH
  cd ..
done
