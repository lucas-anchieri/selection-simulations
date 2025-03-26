#!/bin/bash

# simple script to list paths to all 216 files (each containing 10,000 replicates of a
# single selection scenario) into one file

main_dir=$(pwd)
work_dir=${main_dir}/data/simulation_Ne10000_Ngen1500
cd $work_dir

find "$PWD" -name "*.tsv" > simfiles.txt
