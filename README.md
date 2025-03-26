This code can be used to replicate the analyses in Assessing Ancient DNA Sampling Strategies for Natural Selection Inference in Humans using Allele Frequency Time-Series Data. A preprint of this paper can be found here.

A few disclaimers first. I tried to clean up and comment these scripts as well as I can. There might be a few discrepancies in the file paths that I couldn't catch, but otherwise they should be sufficient to replicate the analyses without too much dificulty. These scripts kind of evolved organically (i.e., messily) over the years, and I had to find some tricks to deal with the limitations of the queuing system of our HPC, so there probably are ways to code some of these scripts that are more optimized and more elegant. Still, I consider that they are worth sharing. Finally, we use four methods to estimate selection in our framework: ApproxWF, BMWS, Slattice, and Sr. In order to run the estimations, one should first install said methods and add them to their `$PATH`.

In the paper, we simulate allele trajectories for a range of selection coefficients and only with starting frequencies of 10% at the first generation. These scripts, however, also allow to perform the simulations for a range of initial frequencies.

I will now go over the different scripts and how they are used.

# Perform the actual simulations
## 01
## 02
## 03
# Sample the simulated trajectories
## 04
## sampling script
## sample files
# Run the methods to estimate selection
## 05
## 06
## 07
## 08
# Plot the results
## 09
