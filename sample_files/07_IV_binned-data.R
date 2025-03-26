#!/usr/bin/env Rscript

# loading libraries
libraries <- c("data.table")
suppressMessages(sapply(libraries, require, character.only = TRUE))

# get current working directory
cwd <- paste0(getwd(), "/")

# set values for initial allele frequency and selection coefficient
initvals <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
              0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
svals <- c(0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)

# define the new time points after binning the data
# generation time
new_gens <- as.character(c((1 + 27) / 2 ,
                           (109 + 111) / 2,
                           round((163 + 181 + 201 + 221 + 232*2 + 235) / 7),
                           round((296 + 301 + 303*2 + 308 + 311*3 + 318*4 + 319*3 + 321*4 + 322*9) / 28),
                           341))

# sample sizes
new_ns <- as.character(c(2, 2, 7, 28, 50))

# for each combination of initial freq. and s
for (init in initvals) {
  for (s in svals) {
    # set file name for sample file with all replicates of sampling
    samplefile <- paste0(wd, "init", init, "/s", s,
                         "/init", init, "_s", s, "_04_II_ancient-like.tsv")

    # set the  file name for the corresponding binned data
    binfile <- paste0(wd, "init", init, "/s", s,
                      "/init", init, "_s", s,
                      "_07_IV_binned-data.tsv")

    # get the sampling data
    samples <- fread(file = samplefile, data.table = FALSE)
    samples[, 1:10] # sanity check

    # create a dataframe for the binned data
    samples_binned <- data.frame(gen = new_gens,
                                 n = new_ns)

    # for each replicates of the sampling data
    for (col in seq(3, dim(samples)[2])) {
      # get the samples in a vector
      d <- samples[, col]

      # bin the samples into fewer time points
      d_binned <- data.frame(d = c(d[1] + d[2],
                                   d[3] + d[4],
                                   d[5] + d[6] + d[7] + d[8] + d[9] + d[10],
                                   d[11] + d[12] + d[13] + d[14] + d[15] + d[16] + d[17] + d[18] + d[19],
                                   d[20]))

      # add to the dataframe
      samples_binned <- cbind(samples_binned, d_binned)
    }

    # we now have a dataframe with as many columns as replicates that we can save
    # with the rest of the sampling data
    column_names <- c("gen", "n", rep("d", dim(samples)[2] - 2))
    write.table(rbind(column_names, samples_binned), file = binfile,
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
}