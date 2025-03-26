#!/usr/bin/env Rscript

# given a file containing simulated allele counts and a file containing a sampling scheme,
# this script samples each replicate trajectory accordingly


# parsing through the arguments
args = commandArgs(trailingOnly=TRUE)

params = list(simfile = "",
              nsamp = 10,
              sampsize = 1,
              interval = 1,
              maxrep = 0,
              sampfile = "sampling_scheme.tsv",
              suffix = format(Sys.time(), "%Y%m%e_%H%M"))

current_arg <- "simfile"
for (arg in args) {
  names(params)
  if (substring(arg,1,2) == "--"){
    current_arg <- strsplit(arg,"--")[[1]][2]
    if (current_arg == "nsamp") {
      samp_mode <- "nsamp"
    } else if (current_arg == "interval") {
      samp_mode <- "interval"
    } else if (current_arg == "sampfile") {
      samp_mode <- "file"
    }
  } else {
    params[current_arg] <- arg
  }
}

# load the simulated allele trajectories
simreps <- read.table(params$simfile, header = F)

# get number of generations and number of replicates from the dimensions
ngen <- dim(simreps)[1]
nrep <- dim(simreps)[2]

# get population size from info stored in the path
path_sub <- strsplit(params$simfile, "_Ne")[[1]][2]
TwoNe <- as.numeric(strsplit(path_sub, "_Ngen")[[1]][1]) * 2

# in our case, the sampling mode is always from a file
if (samp_mode == "file") {
  # load sampling scheme
  sampled_gens <- read.table(params$sampfile, header = TRUE)
} else if (samp_mode == "interval"){
  # this one allows to sample at regular intervals without having to specify a file
  interval <- as.numeric(params$interval)
} else if (samp_mode == "nsamp") {
  # this one allows to keep the same sample size at every gen. without needing a file
  interval <- ngen/as.numeric(params$nsamp)
}

# if no file is provided, create dataframe from scratch
if (samp_mode != "file"){
  sampled_gens <- data.frame(gen = numeric(), n = numeric())
  gen <- interval
  while (gen <= ngen) {
    sample <- list(gen = gen, n = as.numeric(params$sampsize))
    sampled_gens <- rbind(sampled_gens, sample)
    gen <- gen + interval
  }
  sampled_gens$gen <- round(sampled_gens$gen)
}

# select replicates where the derived allele has not been lost
surviving <- which(as.numeric(simreps[ngen,]) != 0)

samples <- sampled_gens

if (params$maxrep == 0){
  params$maxrep <- nrep
}
rep <- 1
for (i in (surviving)) { # for each replicate with non-lost allele
  d <- c()
  for (gen in sampled_gens$gen) { # for each sampled generation the the sampling scheme
    # drawing n individuals at random (specified in file) with allele freq. as probability
    d <- c(d, rbinom(1, sampled_gens[which(sampled_gens$gen==gen),"n"],
                     simreps[gen,i]/TwoNe))
  }
  
  # store the replicates in one object
  samples <- cbind(samples, d)
  
  # if the desired amount of replicates has been reached, stop
  rep <- rep + 1
  if (rep > as.numeric(params$maxrep)){
    break
  }
}

# set name of output file from intial freq. and s information
file_prefix <- strsplit(params$simfile, "/")[[1]]
init <- as.numeric(strsplit(strsplit(path_sub, "init")[[1]][2], "/")[[1]][1])
sval <- as.numeric(
  strsplit(strsplit(strsplit(path_sub,
                             ".tsv")[[1]][1],
                    "_s")[[1]][3],
           "_")[[1]][1])
out_prefix <- paste0("init", init, "_s", sval, "_", params$suffix)
outname <- paste0(out_prefix, ".tsv")

# save the sampling data
write.table(samples, outname, row.names=FALSE, sep="\t", quote = FALSE)