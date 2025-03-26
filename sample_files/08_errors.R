#!/usr/bin/env Rscript

# loading libraries
libraries <- c("data.table")
suppressMessages(sapply(libraries, require, character.only = TRUE))

# get current working directory
cwd <- paste0(getwd(), "/")

# set values for initial allele frequency and selection coefficient
# initvals <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
#               0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
# svals <- c(0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
initvals <- 0.1
svals <- 0.02

# set values for error rates
rate_vals <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

# for each combination of initial freq. and s
for (init in initvals) {
  for (s in svals) {
    # set file name for sample file with all replicates of sampling
    samplefile <- paste0(wd, "init", init, "/s", s,
                         "/init", init, "_s", s, "_04_II_ancient-like.tsv")
    
    # get the sampling data
    samples <- fread(file = samplefile, data.table = FALSE)
    
    # for each error rate
    for (rate in rate_vals) {
      # set the  file name for the corresponding data with errors
      damfile <- paste0(wd, "init", init, "/s", s,
                        "/init", init, "_s", s,
                        "_08_damage_", rate, ".tsv")
      
      # create a dataframe for the data with errors
      samples_damage <- data.frame(gen = samples$gen,
                                   n = samples$n)
      
      # for each replicates of the sampling data
      for (col in seq(3, dim(samples)[2])) {
        # get the samples in a smaller dataframe
        sample <- samples[, c(1, 2, col)]
        
        # set vector to store samples after adding errors
        d_damage <- c()
        
        # for each time point in the current replicate sampling data
        for (row in seq(dim(sample)[1])) {
          # get allele counts
          n <- sample[row, "n"]
          d <- sample[row, "d"]
          
          # set the new derived allele count
          d_new <- d
          if (d != 0) {
            # if there are derived alleles, remove some at random with
            # a probability equal to the error rate
            d_new <- d_new - rbinom(1, d, rate)
          }
          
          if (n-d != 0) {
            # if there are ancestral alleles, remove some at random (aka add derived) with
            # a probability equal to the error rate
            d_new <- d_new + rbinom(1, n-d, rate)
          }
          
          # add the new allele count to the overall vector for the current replicate
          d_damage <- c(d_damage, d_new)
        }
        
        # gather samples with errors from all replicates into a dataframe
        d_damage.df <- data.frame(d = d_damage)
        samples_damage <- cbind(samples_damage, d_damage.df)
      }
      
      # we now have a dataframe with as many columns as replicates that we can save
      # with the rest of the sampling data
      column_names <- c("gen", "n", rep("d", dim(samples)[2] - 2))
      write.table(rbind(column_names, samples_damage), file = damfile,
                  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
}


