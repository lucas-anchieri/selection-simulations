#!/bin/bash -l

# this is a wrapper script that launches 05b_estimate.sh for each combination of selection
# coefficient and initial frequency. Again, due to queuing management limits, we had to
# run it several times in batches. The sample variable has to be changed manually and
# the present script rerun for each dataset/sampling scheme

#s=0.02
#init=0.1
sample="01_Ia_ideal_100g"

#s values: 0 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
#init values: 0.01 0.05 0.1 0.25 0.5

for init in 0.01 0.05 0.1 0.25 0.5
do
	for s in 0 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
	do
		sbatch 05b_estimate.sh $init $s $sample
	done
done