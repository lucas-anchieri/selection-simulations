#!/usr/bin/env Rscript

# loading libraries
libraries <- c("ggplot2", "RColorBrewer", "patchwork")
suppressMessages(sapply(libraries, require, character.only = TRUE))

# get current working directory
cwd <- paste0(getwd(), "/")

# function to isolate ggplot legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#--------------------------------------------------------#
#          Plot the Ideal Datasets (Ia, Ib, Ic)          #
#--------------------------------------------------------#

# list of dataset files
files <- c("01_Ia_ideal_100g.tsv",
           "02_Ib_ideal_500g.tsv",
           "03_Ic_ideal_1000g.tsv",
           "04_II_ancient-like.tsv",
           "05_IIIa_even-timepoints.tsv",
           "06_IIIb_even-samplesizes.tsv")

# starting with the ideal datasets (Ia, Ib, Ic)
# load the sampling schemes
dataset_01 <- read.table(paste0(cwd, "../sample_files/", files[1]),
                         header = TRUE, sep = "\t")
dataset_01$label <- "Dataset Ia - 100 Gen."

dataset_02 <- read.table(paste0(cwd, "../sample_files/", files[2]),
                         header = TRUE, sep = "\t")
dataset_02$label <- "Dataset Ib - 500 Gen."

dataset_03 <- read.table(paste0(cwd, "../sample_files/", files[3]),
                         header = TRUE, sep = "\t")
dataset_03$label <- "Dataset Ic - 1000 Gen."

dataset_I <- rbind(dataset_01, dataset_02, dataset_03)

# compute cumulative sample size from allele counts
samples_cumul <- dataset_I
samples_cumul$count <- 0

for (label in unique(dataset_I$label)) {
  rows <- which(dataset_I$label == label)
  
  count <- 0
  for (row in rows) {
    generation <- dataset_I[row, "gen"]
    count <- count + dataset_I[row, "n"]
    samples_cumul[row, "count"] <- count
  }
}

# choose nice colors
cols <- brewer.pal(12,"Set3")[c(4, 5, 8)]

# plot
IdealDatasets_tmp <- ggplot(data = samples_cumul,
                            aes(x=gen, y=count, fill = label, shape = label)) +
  geom_line(aes(col = label), linewidth = 0.5, alpha=1, linetype = "solid") +
  geom_point(size=4, stroke=0.25, alpha=1, col = "white", fill = "white") +
  geom_point(size=3, stroke=0.25, alpha=1, col = "black") +
  scale_x_continuous(minor_breaks = seq(0, 1000, 50),
                     breaks = seq(0, 1000, 250),
                     name = "Generations",
                     limits = c(0, 1010)) +
  scale_y_continuous(minor_breaks = seq(0, 1060, 100),
                     breaks = seq(0, 1060, 500),
                     name = "Cum. Sample Size",
                     limits = c(0, 1010)) +
  scale_fill_manual(values = cols, name = "Ideal Datasets") +
  scale_color_manual(values = cols, name = "Ideal Datasets") +
  scale_shape_manual(values = c(21, 22, 24), name = "Ideal Datasets") +
  ggtitle("A") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.title.position = "plot",
    plot.title = element_text(margin = margin(b = -15, l = -5)),
    panel.grid = element_line(color = "grey88"),
    legend.text = element_text(size=10),
    legend.position = "right", 
    legend.direction = "vertical",
    legend.title = element_text(size=10, vjust = .5, hjust = 0),
    legend.box = "vertical",
    legend.justification = "left",
    legend.box.just = "left",
    plot.background = element_rect(fill='transparent', color=NA),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white", color=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color="black", linewidth = 0.25),
    axis.title = element_text(size = 10)
  )

IdealDatasets_tmp

# remove the legend and store it on the side
legend <- g_legend(IdealDatasets_tmp)
IdealDatasets_tmp <- IdealDatasets_tmp + theme(legend.position = "none")

# add the legend again with patchwork for better control of the layout
IdealDatasets <- IdealDatasets_tmp + legend + plot_layout(widths = c(7, 3))
IdealDatasets

# save the plot
ggsave(plot = IdealDatasets,
       filename = paste0(cwd, "../plots/panels/Figure_01/01-samples_a_IdealDatasets.pdf"),
       width = 4.5, h = 2, bg='transparent',
       scale = 1.5)

# keep the legend
legend_ideal <- legend

#------------------------------------------------------------------------#
#          Plot the Ancient-like (II), Regularity (IIIa, IIIb),          #
#                        and Binned (IV) Datasets                        #
#------------------------------------------------------------------------#

# load the files
dataset_04 <- read.table(paste0(cwd, "../sample_files/", files[4]),
                         header = TRUE, sep = "\t")
dataset_04$label <- "Dataset II - Ancient-like Data"

dataset_05 <- read.table(paste0(cwd, "../sample_files/", files[5]),
                         header = TRUE, sep = "\t")
dataset_05$label <- "Dataset IIIa - Equal Timing"

dataset_06 <- read.table(paste0(cwd, "../sample_files/", files[6]),
                         header = TRUE, sep = "\t")
dataset_06$label <- "Dataset IIIb - Equal Sample Sizes"

# recreate the binned sampling scheme
dataset_07 <- data.frame(gen = c(14, 110, 206, 311, 341),
                         n = c(2, 2, 7, 28, 50))
dataset_07$label <- "Dataset IV - Binned Data"

datasets_II_IV <- rbind(dataset_04, dataset_05, dataset_06, dataset_07)

datasets_II_IV$label <- factor(datasets_II_IV$label,
                               levels = c("Dataset II - Ancient-like Data",
                                          "Dataset IIIa - Equal Timing",
                                          "Dataset IIIb - Equal Sample Sizes",
                                          "Dataset IV - Binned Data"))

# compute cumulative sample size from allele counts
samples_cumul <- datasets_II_IV
samples_cumul$count <- 0

for (label in unique(datasets_II_IV$label)) {
  rows <- which(datasets_II_IV$label == label)
  
  count <- 0
  for (row in rows) {
    generation <- datasets_II_IV[row, "gen"]
    count <- count + datasets_II_IV[row, "n"]
    samples_cumul[row, "count"] <- count
  }
}

# choose another set of nice colors
cols <- brewer.pal(12,"Set3")[c(4, 5, 8, 11)]

# set plot limits and guides breaks
lims <- c(0,341)
mbrks <- seq(lims[1] + 11, lims[2] + 1, 10)
brks <- c(0, seq(lims[1] + 41, lims[2] + 1, 50))
lbls <- paste0(brks - 341, "  ")
lbls[8] <- "0"
lbls

# plot
AncientDataset_tmp <- ggplot(data = samples_cumul,
                             aes(x=gen, y=count, fill = label, shape = label)) +
  geom_line(aes(col = label), linewidth=0.5, alpha=1, linetype = "solid") +
  geom_point(size=4, stroke=0.25, alpha=1, col = "white", fill = "white") +
  geom_point(size=3, stroke=0.25, alpha=1, col = "black") +
  scale_x_continuous(minor_breaks = mbrks,
                     breaks = brks,
                     labels = lbls,
                     name = "Generations") +
  scale_y_continuous(minor_breaks = seq(0, 100, 10),
                     breaks = seq(0, 100, 40),
                     name = "Cum. Sample Size") +
  scale_fill_manual(values = cols, name = "Imperfect Datasets") +
  scale_color_manual(values = cols, name = "Imperfect Datasets") +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Imperfect Datasets") +
  ggtitle("B") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_line(color = "grey88"),
    plot.title.position = "plot",
    plot.title = element_text(margin = margin(b = -15, l = -5)),
    legend.text = element_text(size=10),
    legend.position = "right", 
    legend.direction = "vertical",
    legend.title = element_text(size=10, vjust = .5, hjust = 0),
    legend.box = "vertical",
    legend.justification = "left",
    legend.box.just = "left",
    plot.background = element_rect(fill='transparent', color=NA),
    legend.box.background = element_rect(fill = "white", color=NA),
    legend.key = element_rect(fill = "white", color=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color="black", linewidth = 0.25),
    axis.title = element_text(size = 10)
  )

AncientDataset_tmp

# remove legend
legend <- g_legend(AncientDataset_tmp)
AncientDataset_tmp <- AncientDataset_tmp + theme(legend.position = "none")

# add the legend again with patchwork for better control of the layout
AncientDataset <- AncientDataset_tmp + legend +
  plot_layout(widths = c(6, 4))
AncientDataset

# save the plot
ggsave(plot = AncientDataset,
       filename = paste0(cwd, "../plots/panels/Figure_01/01-samples_b_AncientDatasets_all.pdf"),
       width = 4.5, h = 2, bg='transparent',
       scale = 1.5)

# keep the legend
legend_ancient <- legend

#-------------------------------------------#
#          Plot the final Figure 1          #
#-------------------------------------------#

# plotting both together with patchwork and saving the final Figure 1
samples <- (IdealDatasets_tmp + legend_ideal +
  AncientDataset_tmp + legend_ancient) +
  plot_layout(ncol = 2, nrow = 2, widths = c(6, 4))
samples

ggsave(plot = samples,
       filename = paste0(cwd, "../plots/Figure_01_sampling.pdf"),
       width = 6, h = 3.5, bg='transparent',
       scale = 1.2)
