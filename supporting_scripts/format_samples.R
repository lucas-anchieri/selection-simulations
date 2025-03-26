#!/usr/bin/env Rscript

# takes simple file containing sampling data and outputs it in several files that are
# each compatible with a different method

# get file path
sampleFile <- commandArgs(trailingOnly=TRUE)[1]

# load it
sample <- read.table(sampleFile, header = TRUE)

# set output file prefix
out_prefix <- "sample_"

#### APPROXWF ####
# set output file name
out_approxwf <- paste0(out_prefix, "approxwf.loci")

# set a single line with name L1 and sampling data for each time point in the format
# "[derived allele count]/[total allele count]"
L1 = c("L1", paste(sample$d, sample$n, sep = "/"))
sample_approxwf <- t(matrix(L1))

# set the column names as dates of the time points in generations
colnames(sample_approxwf) <- c("time", sample$gen)

# write it out
write.table(sample_approxwf, file = out_approxwf, sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE,
            append = FALSE)

#### SLATTICE ####
# set output file name
out_slattice <- paste0(out_prefix, "slattice.txt")

# transposing the matrix is enough for this one
sample_slattice <- t(sample)

# write it out
write.table(sample_slattice, file = out_slattice, sep = ",",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

#### SR ####
# set output file name
out_sr <- paste0(out_prefix, "sr.txt")

# translate time point dates from generations to years backwards in time
timeyears <- (sample$gen - max(sample$gen))*25

# dispose the data in a dataframe
sample_sr <- data.frame(d = sample$d, n = sample$n,
                        t1 = timeyears, t2 = timeyears)

# write it out
write.table(sample_sr, file = out_sr, sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

#### BMWS ####
# for this one, we have to convert the allele count data to .vcf and .meta files

# set output files names
out_bmws <- paste0(out_prefix, "bmws.vcf")
out_bmws_simple <- paste0(out_prefix, "bmws_simple.txt")

# since it needs a .vcf, create a placeholder start of line. the data here will not
# actually be used
genotypes <- c("1", "1", "var", "G", "A", "100", "PASS", ".", "GT")

# create a dataframe for the simple output
sample_bmws_simple <- data.frame()
all_gen <- 1

# create a dataframe for the individuals metadata
meta <- data.frame()

# for each line in the sample file
for (i in seq(dim(sample)[1])) {
  # get generation, date in years, allele counts, genotypes
  gen <- sample[i, "gen"]
  date <- 25 * (max(sample$gen) - gen)
  n <- sample[i, "n"]
  d <- sample[i, "d"]
  genotypes <- c(genotypes, rep("1/1", d))
  genotypes <- c(genotypes, rep("0/0", n-d))
  
  # create placeholder metadata for all individuals of the current time point
  for (ind in seq(n)) {
    indID <- paste(i, ind, sep = "_") # here, we set IDs as "[gen]_[i]"
    ind_data <- data.frame(ID = indID,
                           DateBP = date,
                           Region = "SomePlace",
                           Latitude = 10,
                           Longitude = 10)
    meta <- rbind(meta, ind_data)
  }
  
  # creating a "simple" output that might be necessary later on
  # this one requires one line per generation, even those that were not sampled
  
  # bridge gap between last time point and current time point with 0 at each gen
  while (all_gen != gen) {
    sample_bmws_simple <- rbind(sample_bmws_simple,
                                data.frame(gen = all_gen, n = 0, d = 0))
    all_gen <- all_gen + 1
  }
  
  # for the current time point, add the actual data
  sample_bmws_simple <- rbind(sample_bmws_simple,
                              data.frame(gen = all_gen, n = n, d = d))
  all_gen <- all_gen + 1
}

# output the individuals metadata
colnames(meta)[1] = "#ID"
write.table(meta, file = "sample_bmws.meta", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# create the end of the first line in the .vcf with individuals IDs
vcf_line <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
              "FORMAT", meta$`#ID`)

# write it to the already existing empty .vcf file
write.table(t(vcf_line), file = out_bmws, sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE,
            append = TRUE)

# append the genotype data into the .vcf
sample_bmws <- t(matrix(genotypes))
write.table(sample_bmws, file = out_bmws, sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE,
            append = TRUE)

# also save the "simple" formatting the we produced with one line per generation
write.table(sample_bmws_simple[, -1], file = out_bmws_simple, sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE,
            append = TRUE)
