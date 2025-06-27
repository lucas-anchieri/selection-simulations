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

#-----------------------------------------------------#
#          Load all replicates into "repres"          #
#-----------------------------------------------------#

# loading base ancient-like results
repres_04 <- read.table(paste0(cwd, "../results/04_II_ancient-like_allrepres.tsv"),
                        header = TRUE)

to_keep <- which((repres_04$init == 0.1) & (repres_04$s_true == 0.02))
repres_04 <- repres_04[to_keep, ]
repres_04$label <- "Ancient-like data"

# loading ancient-like results with regular time points
repres_05 <- read.table(paste0(cwd, "../results/05_IIIa_even-timepoints_allrepres.tsv"),
                        header = TRUE)

to_keep <- which((repres_05$init == 0.1) & (repres_05$s_true == 0.02))
repres_05 <- repres_05[to_keep, ]
repres_05$label <- "Even Timing"

# loading ancient-like results with even sample sizes
repres_06 <- read.table(paste0(cwd, "../results/06_IIIb_even-samplesizes_allrepres.tsv"),
                        header = TRUE)

to_keep <- which((repres_06$init == 0.1) & (repres_06$s_true == 0.02))
repres_06 <- repres_06[to_keep, ]
repres_06$label <- "Even Sample Sizes"

# loading binned ancient-like results
repres_07 <- read.table(paste0(cwd, "../results/07_IV_binned-data_allrepres.tsv"),
                        header = TRUE)

to_keep <- which((repres_07$init == 0.1) & (repres_07$s_true == 0.02))
repres_07 <- repres_07[to_keep, ]
repres_07$label <- "Binned Data"

# merging
repres <- rbind(repres_04, repres_05, repres_06, repres_07)

# do some cleaning up
repres[repres$method == "approxwf", "method"] <- "ApproxWF"
repres[repres$method == "bmws", "method"] <- "BMWS"
repres[repres$method == "slattice", "method"] <- "Slattice"
repres[repres$method == "sr", "method"] <- "Sr"

repres$method <- factor(repres$method, levels = c("ApproxWF",
                                                  "BMWS",
                                                  "Slattice",
                                                  "Sr"))

repres$label <- factor(repres$label, levels = c("Ancient-like data",
                                                "Even Timing",
                                                "Even Sample Sizes",
                                                "Binned Data"))

#-------------------------------------------------#
#          Load summary stats into "res"          #
#-------------------------------------------------#

# loading base ancient-like results
res_04 <- read.table(paste0(cwd, "../results/04_II_ancient-like_allres.tsv"),
                     header = TRUE)

to_keep <- which((res_04$init == 0.1) & (res_04$s_true == 0.02))
res_04 <- res_04[to_keep, ]
res_04$label <- "Ancient-like data"

# loading ancient-like results with regular time points
res_05 <- read.table(paste0(cwd, "../results/05_IIIa_even-timepoints_allres.tsv"),
                     header = TRUE)

to_keep <- which((res_05$init == 0.1) & (res_05$s_true == 0.02))
res_05 <- res_05[to_keep, ]
res_05$label <- "Even Timing"

# loading ancient-like results with even sample sizes
res_06 <- read.table(paste0(cwd, "../results/06_IIIb_even-samplesizes_allres.tsv"),
                     header = TRUE)

to_keep <- which((res_06$init == 0.1) & (res_06$s_true == 0.02))
res_06 <- res_06[to_keep, ]
res_06$label <- "Even Sample Sizes"

# loading binned ancient-like results
res_07 <- read.table(paste0(cwd, "../results/07_IV_binned-data_allres.tsv"),
                     header = TRUE)

to_keep <- which((res_07$init == 0.1) & (res_07$s_true == 0.02))
res_07 <- res_07[to_keep, ]
res_07$label <- "Binned Data"

# merging
res <- rbind(res_04, res_05, res_06, res_07)

# do some cleaning up
res$label <- factor(res$label, levels = c("Ancient-like data",
                                          "Even Timing",
                                          "Even Sample Sizes",
                                          "Binned Data"))

res[res$method == "approxwf", "method"] <- "ApproxWF"
res[res$method == "bmws", "method"] <- "BMWS"
res[res$method == "slattice", "method"] <- "Slattice"
res[res$method == "sr", "method"] <- "Sr"

res$method <- factor(res$method, levels = c("ApproxWF",
                                            "BMWS",
                                            "Slattice",
                                            "Sr"))

#----------------------------------------------#
#          Computing additional stats          #
#----------------------------------------------#

repres_stats <- data.frame()
head(repres)
for (method in unique(repres$method)) {
  for (s_true in unique(repres$s_true)) {
    for (label in unique(repres$label)) {
      
      temp <- repres[which(repres$method == method), ]
      temp <- temp[which(temp$s_true == s_true), ]
      temp <- temp[which(temp$label == label), ]
      
      variance <- var(temp$est, na.rm=TRUE)
      bias <- mean(temp$est - s_true, na.rm=TRUE)
      quantiles <- quantile(temp$est, probs = seq(0, 1, 0.05), na.rm=TRUE)
      quantile50 <- as.numeric(quantiles["75%"] - quantiles["25%"])
      quantile90 <- as.numeric(quantiles["95%"] - quantiles["5%"])
      var90 <- var(sort(temp$est)[seq(1, 0.9*length(temp$est))],
                   na.rm = TRUE)
      biasXvar <- sqrt((variance**2)+(bias**2))
      
      stats <- data.frame(method = method,
                          s_true = s_true,
                          label = label,
                          var = variance,
                          var90 = var90,
                          bias = bias,
                          quantile50 = quantile50,
                          quantile90 = quantile90,
                          biasXvar)
      repres_stats <- rbind(repres_stats, stats)
    }
  }
}

# clean up
repres_stats$method <- factor(repres_stats$method, levels = c("ApproxWF",
                                                              "BMWS",
                                                              "Slattice",
                                                              "Sr"))

repres_stats$label <- factor(repres_stats$label, levels = c("Ancient-like data",
                                                            "Even Timing",
                                                            "Even Sample Sizes",
                                                            "Binned Data"))

head(repres_stats)

#----------------------------------------------#
#          Reducing to only ApproxWF           #
#----------------------------------------------#

repres <- repres[repres$method == "ApproxWF", ]
res <- res[res$method == "ApproxWF", ]
repres_stats <- repres_stats[repres_stats$method == "ApproxWF", ]

#----------------------#
# Plot point estimates #
#----------------------#

# set position *pos* of each method
pos <- 0
repres <- cbind(repres, rep(0,dim(repres)[1]))
colnames(repres)[length(colnames(repres))] <- "pos"
for (meth in sort(unique(repres$method))) {
  repres[which(repres$method == meth),"pos"] <- pos
  pos <- pos + 1
}

# plot and save!!
point_estimates <- ggplot(data = repres,
                          aes(x=label, y=est,
                              color="replicates", fill ="replicates")) +
  facet_grid(, vars(method)) +
  geom_point(shape = 16, alpha=0.7, size= 0.5, position = position_jitterdodge(jitter.width = 0.7)) +
  geom_boxplot(alpha = 0, col = "grey25", outlier.shape = NA, linewidth = 0.25) +
  geom_segment(aes(x=0.25, xend = 4.75, y = s_true, yend = s_true),
               color="#E41A1C", linewidth = 0.5, linetype = "dashed") + 
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "grey") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Dataset", y = "Estimated s", fill = "", color = "") +
  coord_cartesian() +
  ggtitle("A") +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(color = "grey95"),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size=8),
        legend.direction = "vertical",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "vertical",
        legend.justification = "center",
        legend.box.just = "left",
        plot.background = element_rect(fill='white', color=NA),
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "none",
        legend.key = element_rect(fill = "white", color=NA),
        axis.title = element_text(size = 10),
        axis.line = element_line(color="black", linewidth = 0.25),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
  )

point_estimates

ggsave(point_estimates,
       file = paste0(cwd, "../plots/panels/Figure_07/07-V-VII_a_point_estimates.png"),
       width = 2.5, height = 2.5, dpi = 600, scale = 1.2)

#------------------------#
#          RMSE          #
#------------------------#

rmse <- ggplot(data = res,
               aes(x=label, y=rmse, col = as.factor(method))) +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.01, 0.005),
                     minor_breaks = seq(0, 0.01, 0.001),
                     limits = c(0, 0.01)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,3,2,4)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,3,2,4)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24, 23),
                     name = "Method") +
  xlab("Dataset") +
  ylab("RMSE") +
  ggtitle("C") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size=8),
        legend.direction = "horizontal",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

rmse

ggsave(rmse, file = paste0(cwd, "../plots/panels/Figure_07/07-V-VII_c_rmse.pdf"),
       width = 2.5, height = 2.5, scale = 1.5)

#------------------------#
#          Bias          #
#------------------------#

bias <- ggplot(data = repres_stats,
               aes(x=label, y=bias, col = as.factor(method))) +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(-0.01, 0.0025, 0.001),
                     minor_breaks = seq(-0.01, 0.0025, 0.0005),
                     limits = c(NA, 0.0025)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,3,2,4)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,3,2,4)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24, 23),
                     name = "Method") +
  xlab("Dataset") +
  ylab("Bias") +
  ggtitle("D") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size=8),
        legend.direction = "horizontal",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

bias

ggsave(bias, file = paste0(cwd, "../plots/panels/Figure_07/07-V-VII_d_bias.pdf"),
       width = 2.5, height = 2.5, scale = 1.5)

#----------------------------#
#          Variance          #
#----------------------------#

variance <- ggplot(data = repres_stats,
                   aes(x=label, y=var, col = as.factor(method))) +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.001, 0.00001),
                     minor_breaks = seq(0, 0.001, 0.000005),
                     limits = c(0, NA)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,3,2,4)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,3,2,4)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24, 23),
                     name = "Method") +
  xlab("Dataset") +
  ylab("Variance") +
  ggtitle("E") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(), #element_line(color = "grey95")
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size=8),
        legend.direction = "horizontal",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

variance

ggsave(variance, file = paste0(cwd, "../plots/panels/Figure_07/07-V-VII_e_variance.pdf"),
       width = 2.5, height = 2.5, scale = 1.5)

#-------------------------#
# Plot True Positive Rate #
#-------------------------#

tpr <- ggplot(data = res,
              aes(x=label, y=propPos, col = as.factor(method))) +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 5, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 1, 0.5),
                     minor_breaks = seq(0, 1, 0.1),
                     labels = paste0(seq(0, 100, 50), "%"),
                     limits = c(0, 1)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24),
                     name = "Method") +
  xlab("Dataset") +
  ylab("True Positive Rate") +
  ggtitle("B") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -10, l = -5)),
        legend.text = element_text(size=8),
        legend.direction = "horizontal",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.title = element_text(size = 10),
        axis.title.y = element_text(hjust = 0),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

tpr

ggsave(tpr, file = paste0(cwd, "../plots/panels/Figure_07/07-V-VII_b_tpr.pdf"),
       width = 2.5, height = 2.5, scale = 1.5)

#----------------------------#
#       FINAL FIGURE 6       #
#----------------------------#

# removing legends
legend <- g_legend(rmse)
rmse <- rmse + theme(legend.position = "none")
bias <- bias + theme(legend.position = "none")
tpr <- tpr + theme(legend.position = "none")
variance <- variance + theme(legend.position = "none")

# removing x axis on point-estimates, rmse, and bias
point_estimates <- point_estimates + theme(axis.title.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.text.x = element_blank())
rmse <- rmse + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())
bias <- bias + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())

# approxwf only one figure for all
fig_V.VII_1 <- (((point_estimates / tpr) + plot_layout(heights = c(2, 1))) |
                  ((rmse / bias / variance)) + plot_layout(heights = c(1, 1, 1))) /
  legend + plot_layout(heights = c(9.8, 0.2))

fig_V.VII_1

ggsave(fig_V.VII_1, file = paste0(cwd, "../plots/Figure_07_regularity.pdf"),
       width = 5, height = 5, scale = 1.2)
