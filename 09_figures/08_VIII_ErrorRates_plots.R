#!/usr/bin/env Rscript

# loading libraries
libraries <- c("ggplot2", "RColorBrewer", "ggrepel", "patchwork")
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

# values for true s
s_vals_txt <- c("0", "0.01", "0.02", "0.03", "0.04", "0.05", "0.06",
                "0.07", "0.08", "0.09", "0.1")

#-----------------------------------------------------#
#          Load all replicates into "repres"          #
#-----------------------------------------------------#

# loading base ancient-like results
repres_ancient <- read.table(paste0(cwd, "../results/04_II_ancient-like_allrepres.tsv"),
                             header = TRUE)

to_keep <- which((repres_ancient$init == 0.1) & (repres_ancient$s_true == 0.02))
repres_ancient <- repres_ancient[to_keep, ]
repres_ancient$sample <- 0

# loading results with error rates
repres01 <- read.table(paste0(cwd, "../results/08_damage_0.01_allrepres.tsv"),
                        header = TRUE)
repres01$sample <- 0.01

repres05 <- read.table(paste0(cwd, "../results/08_damage_0.05_allrepres.tsv"),
                       header = TRUE)
repres05$sample <- 0.05

repres10 <- read.table(paste0(cwd, "../results/08_damage_0.1_allrepres.tsv"),
                       header = TRUE)
repres10$sample <- 0.1

repres15 <- read.table(paste0(cwd, "../results/08_damage_0.15_allrepres.tsv"),
                       header = TRUE)
repres15$sample <- 0.15

repres20 <- read.table(paste0(cwd, "../results/08_damage_0.2_allrepres.tsv"),
                       header = TRUE)
repres20$sample <- 0.2

repres25 <- read.table(paste0(cwd, "../results/08_damage_0.25_allrepres.tsv"),
                       header = TRUE)
repres25$sample <- 0.25

repres <- rbind(repres_ancient,
                repres01, repres05, repres10, repres15, repres20, repres25)

# do some cleaning up
repres[repres$method == "approxwf", "method"] <- "ApproxWF"
repres[repres$method == "bmws", "method"] <- "BMWS"
repres[repres$method == "slattice", "method"] <- "Slattice"
repres[repres$method == "sr", "method"] <- "Sr"

repres$method <- factor(repres$method, levels = c("ApproxWF",
                                                  "BMWS",
                                                  "Slattice",
                                                  "Sr"))

#-------------------------------------------------#
#          Load summary stats into "res"          #
#-------------------------------------------------#

# loading base ancient-like results
res_ancient <- read.table(paste0(cwd, "../results/04_II_ancient-like_allres.tsv"),
                          header = TRUE)

to_keep <- which((res_ancient$init == 0.1) & (res_ancient$s_true == 0.02))
res_ancient <- res_ancient[to_keep, ]
res_ancient$sample <- 0

# loading results with errors
res01 <- read.table(paste0(cwd, "../results/08_damage_0.01_allres.tsv"),
                    header = TRUE)
res01$sample <- 0.01

res05 <- read.table(paste0(cwd, "../results/08_damage_0.05_allres.tsv"),
                    header = TRUE)
res05$sample <- 0.05

res10 <- read.table(paste0(cwd, "../results/08_damage_0.1_allres.tsv"),
                    header = TRUE)
res10$sample <- 0.1

res15 <- read.table(paste0(cwd, "../results/08_damage_0.15_allres.tsv"),
                    header = TRUE)
res15$sample <- 0.15

res20 <- read.table(paste0(cwd, "../results/08_damage_0.2_allres.tsv"),
                    header = TRUE)
res20$sample <- 0.2

res25 <- read.table(paste0(cwd, "../results/08_damage_0.25_allres.tsv"),
                    header = TRUE)
res25$sample <- 0.25

res <- rbind(res_ancient,
             res01, res05, res10, res15, res20, res25)

# do some cleaning up
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
    for (samp in unique(repres$sample)) {
      temp <- repres[which(repres$method == method), ]
      temp <- temp[which(temp$s_true == s_true), ]
      temp <- temp[which(temp$samp == samp), ]
      
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
                          sample = samp,
                          var = variance,
                          var90 = var90,
                          bias = bias,
                          quantile50 = quantile50,
                          quantile90 = quantile90,
                          biasXvar = biasXvar)
      repres_stats <- rbind(repres_stats, stats)
    }
  }
}

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

y_labs <- c("0%", "1%", "5%", "10%", "15%", "20%", "25%")

point_estimates <- ggplot(data = repres,
                          aes(x=as.factor(sample), y=est,
                              color="replicates", fill ="replicates")) +
  facet_grid(, vars(method)) +
  geom_point(shape = 16, alpha=0.7, size=0.5, position = position_jitterdodge(jitter.width = 0.7)) +
  geom_boxplot(alpha = 0, col = "grey25", outlier.shape = NA, linewidth = 0.25) +
  geom_segment(aes(x=0.25, xend = 7.75, y = s_true, yend = s_true),
               color="#E41A1C", linewidth = 0.5, linetype = "dashed") +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "grey") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Error rate", y = "Estimated s", fill = "", color = "") +
  scale_y_continuous(breaks = seq(-0.07, 0.1, 0.01)) +
  scale_x_discrete(labels = y_labs) +
  coord_cartesian(ylim = c(NA, 0.05)) +
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
       file = paste0(cwd, "../plots/panels/Figure_08/08-VIII_a_point_estimates.pdf"),
       width = 2.6, height = 2.6, scale = 1.2)

#------------------------#
#          RMSE          #
#------------------------#

rmse <- ggplot(data = res,
               aes(x=sample, y=rmse, col = as.factor(method))) +
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
  scale_y_continuous(breaks = seq(0, 0.015, 0.005),
                     minor_breaks = seq(0, 0.015, 0.001),
                     limits = c(0, 0.015)) +
  scale_x_continuous(breaks = c(0, 0.01, seq(0.05, 0.25, 0.05)),
                     minor_breaks = seq(0, 0.25, 0.05),
                     labels = c(0, 0.01, seq(0.05, 0.25, 0.05))*100,
                     limits = c(0, 0.25)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24),
                     name = "Method") +
  xlab("Error Rate (%)") +
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
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

rmse

ggsave(rmse, file = paste0(cwd, "../plots/panels/Figure_08/08-VIII_c_rmse.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#------------------------#
#          Bias          #
#------------------------#

bias <- ggplot(data = repres_stats,
               aes(x=sample, y=bias, col = as.factor(method))) +
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
  scale_y_continuous(breaks = seq(-0.02, 0.02, 0.005),
                     minor_breaks = seq(-0.02, 0.02, 0.001),
                     limits = c(-0.013, 0.005)) +
  scale_x_continuous(breaks = c(0, 0.01, seq(0.05, 0.25, 0.05)),
                     minor_breaks = seq(0, 0.25, 0.05),
                     labels = c(0, 0.01, seq(0.05, 0.25, 0.05))*100,
                     limits = c(0, 0.25)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24),
                     name = "Method") +
  xlab("Error Rate (%)") +
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
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

bias

ggsave(bias, file = paste0(cwd, "../plots/panels/Figure_08/08-VIII_d_bias.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#----------------------------#
#          Variance          #
#----------------------------#

variance <- ggplot(data = repres_stats,
                   aes(x=sample, y=var, col = as.factor(method))) +
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
  scale_y_continuous(breaks = seq(0, 0.008, 0.00001),
                     minor_breaks = seq(0, 0.008, 0.000005),
                     limits = c(0, 0.00005)) +
  scale_x_continuous(breaks = c(0, 0.01, seq(0.05, 0.25, 0.05)),
                     minor_breaks = seq(0, 0.25, 0.05),
                     labels = c(0, 0.01, seq(0.05, 0.25, 0.05))*100,
                     limits = c(0, 0.25)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24),
                     name = "Method") +
  xlab("Error Rate (%)") +
  ylab("Variance") +
  ggtitle("E") +
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
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

variance

ggsave(variance, file = paste0(cwd, "../plots/panels/Figure_08/08-VIII_e_variance.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#-------------------------#
# Plot True Positive Rate #
#-------------------------#

tpr <- ggplot(data = res,
              aes(x=sample, y=propPos, col = as.factor(method))) +
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
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     minor_breaks = seq(0, 1, 0.05),
                     labels = paste0(seq(0, 100, 25), "%"),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.01, seq(0.05, 0.25, 0.05)),
                     minor_breaks = seq(0, 0.25, 0.05),
                     labels = c(0, 0.01, seq(0.05, 0.25, 0.05))*100,
                     limits = c(0, 0.25)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24),
                     name = "Method") +
  xlab("Error Rate (%)") +
  ylab("Power") +
  ggtitle("B") +
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
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

tpr

ggsave(tpr, file = paste0(cwd, "../plots/panels/Figure_08/08-VIII_b_tpr.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#------------------------#
#       FINAL PLOT       #
#------------------------#

# removing legends
legend <- g_legend(rmse)
rmse <- rmse + theme(legend.position = "none")
bias <- bias + theme(legend.position = "none")
tpr <- tpr + theme(legend.position = "none")
variance <- variance + theme(legend.position = "none")

# removing x axis on rmse and bias
rmse <- rmse + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())
bias <- bias + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())

# approxwf only one figure for all
fig_VIII_1 <- (((point_estimates / tpr) + plot_layout(heights = c(2, 1))) |
                  ((rmse / bias / variance)) + plot_layout(heights = c(1, 1, 1))) /
  legend + plot_layout(heights = c(9.8, 0.2))

fig_VIII_1

ggsave(fig_VIII_1, file = paste0(cwd, "../plots/Figure_08_errors.pdf"),
       width = 5, height = 5, scale = 1.2)
