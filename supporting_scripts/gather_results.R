#!/usr/bin/env Rscript

# this scripts gathers all replicate results in separate files into one single file per
# method

library(data.table)

#### APPROXWF ####
# load paths to result files
awf_resfiles <- read.table("approxwf_files.txt")[,1]

# load the first result file
approxwf_temp <- fread(awf_resfiles[1], header = TRUE,
                       data.table = FALSE)

# create dataframe that will store all results, with index column from the first file
est <- data.frame(index = approxwf_temp$index)

# for each result file, take the results and store the posterior as an additional column
i=0
for (resfile in awf_resfiles) {
  # load file
  approxwf_temp <- fread(resfile, header = TRUE,
                         data.table = FALSE)

  # get posterior
  est_temp <- data.frame(est = approxwf_temp$s_L1)
  names(est_temp) <- paste0("est", i)
  
  # add it to the dataframe
  est <- cbind(est, est_temp)

  # increase iterator
  i <- i + 1
}

# save the dataframe with all results in a single file
fwrite(est, file = "approxwfEst.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")

#### SLATTICE ####
# load paths to result files
slattice_resfiles <- read.table("slattice_files.txt")[,1]

# create dataframe that will store all results
est <- data.frame()

# for each result file, take the results and store it as an additional row
for (resfile in slattice_resfiles) {
  # first, the slattice output has to be processed by another script
  system(paste0("./supporting_scripts/process_slattice.py ", resfile))
  
  # load the result file
  slattice_temp <- fread(paste0(resfile, ".simple"), header=F,
                         data.table = FALSE)
  
  # get the point estimate and the confidence interval
  slattice <- data.frame(est = slattice_temp[1,1],
                         minCI = slattice_temp[2,1],
                         maxCI = slattice_temp[3,1])
  
  # add it to the dataframe
  est <- rbind(est, slattice)
}

# save the dataframe with all results in a single file
fwrite(est, file = "slatticeEst.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")

#### SR ####
# load paths to result files

sr_resfiles <- read.table("sr_files.txt")[,1]

# load the first result file
sr_temp <- fread(sr_resfiles[1], header = TRUE,
                 data.table = FALSE)

# create dataframe that will store all results, with index column from the first file
est <- data.frame(index = sr_temp$gen)

n_iter <- length(sr_temp$gen)

# for each result file, take the results and store the posterior as an additional column
i=0
for (resfile in sr_resfiles) {
  # load file
  sr_temp <- fread(resfile, header = TRUE,
                   data.table = FALSE)

  if (dim(sr_temp)[1] < n_iter) {
    # if the chain didn't run to the end, fill the rest with NAs
    est_temp <- data.frame(est = rep(NA, n_iter))
  } else {
    # get posterior
    est_temp <- data.frame(est = sr_temp$alpha2)
  }
  
  # add it to the dataframe
  names(est_temp) <- paste0("est", i)
  est <- cbind(est, est_temp)

  # increase iterator
  i <- i + 1
}

# save the dataframe with all results in a single file
fwrite(est, file = "srEst.tsv",
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t",
       na = "NA")

# #### BMWS ####
# load paths to result files
bmws_resfiles <- read.table("bmws_files.txt")[,1]

# create dataframe that will store all results
est <- data.frame()

# for each result file, take the results and store it as an additional row
for (resfile in bmws_resfiles) {
  # load the result file
  bmws_temp <- fread(resfile, header=F,
                     data.table = FALSE)
  
  # get the point estimate and the confidence interval
  bmws <- data.frame(est = bmws_temp[1,7],
                     minCI = bmws_temp[1,8],
                     maxCI = bmws_temp[1,9])
  
  # add it to the dataframe
  est <- rbind(est, bmws)
}

# save the dataframe with all results in a single file
fwrite(est, file = "bmwsEst.tsv",
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
