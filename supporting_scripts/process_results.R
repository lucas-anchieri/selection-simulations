#!/usr/bin/env Rscript

# computes summary statistics for a given method and scenario and stores everything
# in a single file per scenario

# loading libraries
libraries <- c("hdrcde", "data.table")
suppressMessages(sapply(libraries, require, character.only = TRUE))

# get current working directory
cwd <- paste0(getwd(), "/")

# parsing through the arguments
args = commandArgs(trailingOnly=TRUE)

params = list(resfile = "",
              init = 0.1,
              s_true = 0.02,
              Ne = 10000,
              ngen = 500,
              plotName = "plot.png",
              method = "")

current_arg <- ""
for (arg in args) {
  names(params)
  if (substring(arg,1,2) == "--"){
    current_arg <- strsplit(arg,"--")[[1]][2]
  } else {
    if(is.na(current_arg)) {quit(status=1)}
    params[current_arg] <- arg
  }
}

# prepare variables
nrep <- 0 # number of replicates
npos <- 0 # number of replicates where 0 not in HDR 95%
npos_den <- 0 # number of replicates where positive density > 95%
ntrue <- 0 # number of reps where 0 not in HDR 95% AND s_true in HDR 95%
modes <- c() # modes of the replicates
init <- as.numeric(params$init) # initial allele frequency
s_true <- as.numeric(params$s_true) # selection coefficient

# load the results
print("Loading results from file")
results <- fread(file = paste0(wd, params$resfile), header = TRUE, data.table = FALSE)

### SLATTICE """
if (params$method == "slattice") {
  print("Processing Slattice")
  
  # reading the Slattice replicate results line by line
  for (line in seq(1,dim(results)[1])) {
    row <- results[line, ]
    
    # update total number of replicates
    nrep <- nrep + 1
    
    pos = TRUE # can we significantly reject neutrality and accept positive selection
    strue = FALSE # is the true value of s included in the 95% interval
    
    if (is.na(row$minCI) || is.na(row$maxCI)) {
      # if there is no conficence interval, no significant sign of positive selection
      pos = FALSE
    } else {
      if (0 >= row$minCI && 0 <= row$maxCI) {
        # if 0 is included in the 95% CI, no significant sign of positive selection
        pos = FALSE
      }
      
      if (s_true >= row$minCI && s_true <= row$maxCI) {
        # if s_true is included in the 95% CI
        strue = TRUE
      }
    }
    
    # since we don't have a posterior, pos_den and pos will be the same
    pos_den <- pos
    
    # update counters
    if (pos) {npos <- npos + 1}
    if (pos_den) {npos_den <- npos_den + 1}
    if (strue) {ntrue <- ntrue + 1}
    
    # get mode (aka the point estimate)
    modes <- c(modes, row$est)
    
    # package everything neatly in a dataframe
    repres <- data.frame(Ne = params$Ne,
                         ngen = params$ngen,
                         s_true = s_true,
                         init = init,
                         method = params$method,
                         pos_HDI = pos,
                         pos_den = pos_den,
                         strue = strue,
                         est = row$est)
    
    # write the summary stats for each replicate in the same file, appending a new line
    # in the end, this "repres" file will have as many lines as replicates (1,000)
    write.table(repres, file = paste0(wd, "repres.tsv"),
                row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE,
                col.names = !file.exists(paste0(wd, "repres.tsv")))
    cat("\r", line)
  }
} else if (params$method == "bmws") {
  print("Processing bmws")
  
  # reading the Slattice replicate results line by line
  for (line in seq(1,dim(results)[1])) {
    row <- results[line, ]
    
    # update total number of replicates
    nrep <- nrep + 1
    
    pos = TRUE # can we significantly reject neutrality and accept positive selection
    strue = FALSE # is the true value of s included in the 95% interval
    
    if (is.na(row$minCI) || is.na(row$maxCI)) {
      # if there is no conficence interval, no significant sign of positive selection
      pos = FALSE
    } else {
      if (0 >= row$minCI && 0 <= row$maxCI) {
        # if 0 is included in the 95% CI, no significant sign of positive selection
        pos = FALSE
      }
      
      if (s_true >= row$minCI && s_true <= row$maxCI) {
        # if s_true is included in the 95% CI
        strue = TRUE
      }
    }
    
    # since we don't have a posterior, pos_den and pos will be the same
    pos_den <- pos
    
    # update counters
    if (pos) {npos <- npos + 1}
    if (pos_den) {npos_den <- npos_den + 1}
    if (strue) {ntrue <- ntrue + 1}
    
    # get mode (aka the point estimate)
    modes <- c(modes, row$est)
    
    # package everything neatly in a dataframe
    repres <- data.frame(Ne = params$Ne,
                         ngen = params$ngen,
                         s_true = s_true,
                         init = init,
                         method = params$method,
                         pos_HDI = pos,
                         pos_den = pos_den,
                         strue = strue,
                         est = row$est)
    
    # write the summary stats for each replicate in the same file, appending a new line
    # in the end, this "repres" file will have as many lines as replicates (1,000)
    write.table(repres, file = paste0(wd, "repres.tsv"),
                row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE,
                col.names = !file.exists(paste0(wd, "repres.tsv")))
    cat("\r", line)
  }
} else {
  # if method is ApproxWF or Sr, then we are dealing with a posterior distribution
  print("Processing posteriors")
  
  # removing the first 10% as burn-in
  posteriors <- results[-(seq(1, dim(results)[1]*0.1)), ]
  
  # loop through the posteriors
  cat("\nLoop through posteriors\n")
  for (col in seq(2, dim(posteriors)[2])) {
    if (params$method == "sr") {
      # if method is Sr, the results are given as 2Ne*s and we need to normalize to 2Ne
      posteriors[, col] <- posteriors[, col] / (2 * params$Ne)
    }
    
    # update total number of replicates
    nrep <- nrep + 1
    
    pos = TRUE # can we significantly reject neutrality and accept positive selection
    strue = FALSE # is the true value of s included in the 95% interval
    
    if (!is.na(posteriors[1, col])) {
      # if the posterior is not empty
      # compute Highest Density Regions (HDRs) for the posterior
      dist_hdr <- hdr(den = density(posteriors[, col]))
      
      # isolate the 95% HDR
      intervals95 <- dist_hdr$hdr[2, ]
    } else {
      # if the posterior is empty, fill with NAs
      intervals95 <- c(NA,NA)
      dist_hdr <- data.frame(mode = NA)
      pos = FALSE
    }
    
    # for each 95% HDR interval, check whether 0 and s_true are in it or not
    for (int in seq(1, length(intervals95), 2)) {
      if (is.na(intervals95[int]) || is.na(intervals95[int + 1])) {
        # skipping NAs
        next
      }
      
      if (0 >= intervals95[int] && 0 <= intervals95[int + 1]) {
        # if 0 is included in the 95% HDR, no significant sign of positive selection
        pos = FALSE
      }
      
      if (s_true >= intervals95[int] && s_true <= intervals95[int + 1]) {
        # if s_true is included in the 95% HDR
        strue = TRUE
      }
    }
    
    # check >0 density, aka proportion of the values in the posterior that are above 0
    pos_den <- FALSE 
    interval <- which(posteriors[, col] > 0)
    ZeroDen <- length(interval)/length(posteriors[, col])
    
    if (ZeroDen > 0.95) {
      # if more than 95% of the values in the posterior are above 0, then reject neut.
      # this can be used as an alternative to the CI/HDR to test for significance, but
      # we don't actually use it in the paper
      pos_den <- TRUE
    }
    
    # update counters
    if (pos) {npos <- npos + 1}
    if (pos_den) {npos_den <- npos_den + 1}
    if (strue) {ntrue <- ntrue + 1}
    
    # get mode of the posterior distribution
    modes <- c(modes, dist_hdr$mode)
    
    # package everything neatly in a dataframe
    repres <- data.frame(Ne = params$Ne,
                         ngen = params$ngen,
                         s_true = s_true,
                         init = init,
                         method = params$method,
                         pos_HDI = pos,
                         pos_den = pos_den,
                         strue = strue,
                         est = dist_hdr$mode)
    
    # write the summary stats for each replicate in the same file, appending a new line
    # in the end, this "repres" file will have as many lines as replicates (1,000)
    write.table(repres, file = paste0(wd, "repres.tsv"),
                row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE,
                col.names = !file.exists(paste0(wd, "repres.tsv")))
    cat("\r", col - 1)
  }
}
cat("\nDONE!\n")

# compute RMSE across all modes
print("Computing rmse")
rmse <- sqrt(mean((modes[which(!is.na(modes))] - s_true)^2))

# store summary values across all replicates in a dataframe
print("Storing Result")
res <- data.frame(method = params$method,
                  Ne = params$Ne,
                  init = init,
                  s_true = s_true,
                  nrep = nrep,
                  npos = npos,
                  npos_den = npos_den,
                  ntrue = ntrue,
                  propPos = npos / nrep,
                  propPos_den = npos_den / nrep,
                  propTrue = ntrue / nrep,
                  rmse = rmse
                  )

# write the summary stats over all in the same file, appending a new line
# in the end, this "res" file will have only one line summarizing all replicates
write.table(res, file = paste0(wd, "res.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE,
            col.names = !file.exists(paste0(wd, "res.tsv")))
