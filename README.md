This code can be used to replicate the analyses in Assessing Ancient DNA Sampling Strategies for Natural Selection Inference in Humans using Allele Frequency Time-Series Data. A preprint of this paper can be found here.

A few disclaimers first. I tried to clean up and comment these scripts as well as I can. There might be a few discrepancies in the file paths that I couldn't catch, but otherwise they should be sufficient to replicate the analyses without too much dificulty. These scripts kind of evolved organically (i.e., messily) over the years, and I had to find some tricks to deal with the limitations of the queuing system of our HPC, so there probably are ways to code some of these scripts that are more optimized and more elegant. Still, I consider that they are worth sharing. Finally, we use four methods to estimate selection in our framework: *ApproxWF*, *BMWS*, *Slattice*, and *Sr*. In order to run the estimations, one should first install said methods and add them to their `$PATH`.

In the paper, we simulate allele trajectories for a range of selection coefficients and only with starting frequencies of 10% at the first generation. These scripts, however, also allow to perform the simulations for a range of initial frequencies.

I will now go over the different scripts and how they are used.

# Perform the actual simulations

The first 3 scripts are used to simulate allele trajectories taking into account different parameters: Population size, selection coefficient, initial allele frequency, total number of generations.

## 01_simulate.sh

This script is used to simulate allele trajectories for a range of values for s (from 0 to 0.1) and a range of initial allele frequencies (from 1% to 50%). Since the queueing system of our HPC only allows running 10,000 jobs at a time and we are simulating a lot of replicates, we performed the simulations in 4 batches, by running the script 4 times with different values for s and initial frequency. For each combination of parameters, 10,000 replicates of the allele trajectories are simulated.

## 02_gather_reps.sh

This script gathers all 10,000 replicate of a simulated scenario into a single file.

## 03_list_simfiles.sh

This is a simple script that lists the paths to all 216 files (each containing 10,000 replicates of a single selection scenario) into one file.

# Sample the simulated trajectories

Then, we draw samples at random from the simulated trajectories in order to replicate time-series datasets of allele frequencies.

## 04_sample.sh

This script is a wrapper that launches `supporting_scripts/sampling.R` in parallel for each selection scenario. One of the sampling schemes in the `sample_files/` directory has to be selected, and the script has to be launched separately for each different sampling scheme.

Given a file containing simulated allele counts and a file containing a sampling scheme, `supporting_scripts/sampling.R` will then sample each replicate trajectory accordingly and save them.

# Run the methods to estimate selection

## 05a_launch_estimate.sh and 05b_estimate.sh

`05b_estimate.sh` launches the four methods in order to estimate selection on the time-series samples that were drawn previously. It does so for a given combination of the selection coefficien and initial frequency. `05a_launch_estimate.sh` is a wrapper script that runs `05b_estimate.sh` for each combination of selection coefficient and initial frequency. Those scripts have to be run again for each sampling scheme.

## 06_gather_res.sh

This script launches `supporting_scripts/gather_results.R` in order to list all the result files for all replicates and all methods for a given scenario.

## 07_process_results.sh

This script launches `supporting_scripts/process_results.R`, which computes summary statistics for each method in each scenario, as well as mode and credible interval for ABC methods, and stores everything in a single file per scenario

## 08_gather_all.sh

This script gathers all the result files generated at the last step into two files for a given scenario. One, ending with `all_repres.tsv`, contains independent results for each of the 1,000 replicates for every combination of parameters. The other, ending with `all_res.tsv`, contains summary statistics computed from the replicates for every combination of parameters.

# Plot the results

## 09_figures

The scripts in this folder allow to take the results computed in the last step (and available in the `results/` folder) and plot them into the figures seen in the paper.
