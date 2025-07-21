#!/usr/bin/env Rscript

# loading libraries
libraries <- c("ggplot2", "RColorBrewer", "patchwork", "data.table")
suppressMessages(sapply(libraries, require, character.only = TRUE))

# get current working directory
cwd <- paste0(getwd(), "/")

# function to isolate ggplot legend
g_legend <- function(a.gplot) {
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

repres100 <- read.table(paste0(cwd, "../results/01_Ia_ideal_100g_allrepres.tsv"),
                        header = TRUE)
repres100$sample <- "100 Gen."

repres500 <- read.table(paste0(cwd, "../results/02_Ib_ideal_500g_allrepres.tsv"),
                        header = TRUE)
repres500$sample <- "500 Gen."

repres1000 <- read.table(paste0(cwd, "../results/03_Ic_ideal_1000g_allrepres.tsv"),
                         header = TRUE)
repres1000$sample <- "1000 Gen."

repres <- rbind(repres100, repres500, repres1000)

# do some cleaning up
repres[repres$method == "approxwf", "method"] <- "ApproxWF"
repres[repres$method == "bmws", "method"] <- "BMWS"
repres[repres$method == "slattice", "method"] <- "Slattice"
repres[repres$method == "sr", "method"] <- "Sr"

repres$method <- factor(repres$method, levels = c("ApproxWF",
                                                  "BMWS",
                                                  "Slattice",
                                                  "Sr"))

repres$sample <- factor(repres$sample, levels = c("100 Gen.",
                                                  "500 Gen.",
                                                  "1000 Gen."))

#-------------------------------------------------#
#          Load summary stats into "res"          #
#-------------------------------------------------#

res100 <- read.table(paste0(cwd, "../results/01_Ia_ideal_100g_allres.tsv"),
                     header = TRUE)
res100$sample <- "100 Gen."
res100$sample_old <- "100 Generations"

res500 <- read.table(paste0(cwd, "../results/02_Ib_ideal_500g_allres.tsv"),
                     header = TRUE)
res500$sample <- "500 Gen."
res500$sample_old <- "500 Generations"

res1000 <- read.table(paste0(cwd, "../results/03_Ic_ideal_1000g_allres.tsv"),
                      header = TRUE)
res1000$sample <- "1000 Gen."
res1000$sample_old <- "1000 Generations"

res <- rbind(res100, res500, res1000)

# do some cleaning up
res[res$method == "approxwf", "method"] <- "ApproxWF"
res[res$method == "bmws", "method"] <- "BMWS"
res[res$method == "slattice", "method"] <- "Slattice"
res[res$method == "sr", "method"] <- "Sr"

res$method <- factor(res$method, levels = c("ApproxWF",
                                            "BMWS",
                                            "Slattice",
                                            "Sr"))

res$sample <- factor(res$sample, levels = c("1000 Gen.",
                                            "500 Gen.",
                                            "100 Gen."))

res$sample_rev <- factor(res$sample, levels = c("100 Gen.",
                                                "500 Gen.",
                                                "1000 Gen."))

res$sample_old <- factor(res$sample_old, levels = c("100 Generations",
                                                    "500 Generations",
                                                    "1000 Generations"))

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

      variance <- var(temp$est, na.rm = TRUE)
      bias <- mean(temp$est - s_true, na.rm = TRUE)
      quantiles <- quantile(temp$est, probs = seq(0, 1, 0.05), na.rm = TRUE)
      quantile50 <- as.numeric(quantiles["75%"] - quantiles["25%"])
      quantile90 <- as.numeric(quantiles["95%"] - quantiles["5%"])
      var90 <- var(sort(temp$est)[seq(1, 0.9 * length(temp$est))],
                   na.rm = TRUE)
      biasXvar <- sqrt((variance**2) + (bias ** 2))

      res_temp <- res[which(res$method == method), ]
      res_temp <- res_temp[which(res_temp$s_true == s_true), ]
      res_temp <- res_temp[which(res_temp$sample == samp), ]
      rmse <- res_temp$rmse

      stats <- data.frame(method = method,
                          s_true = s_true,
                          sample = samp,
                          RMSE = rmse,
                          Variance = variance,
                          var90 = var90,
                          Bias = bias,
                          quantile50 = quantile50,
                          quantile90 = quantile90,
                          biasXvar = biasXvar)
      repres_stats <- rbind(repres_stats, stats)
    }
  }
}

# clean up
repres_stats$method <- factor(repres_stats$method, levels = c("ApproxWF",
                                                              "BMWS",
                                                              "Slattice",
                                                              "Sr"))

repres_stats$sample <- factor(repres_stats$sample, levels = c("100 Gen.",
                                                              "500 Gen.",
                                                              "1000 Gen."))

head(repres_stats)

#----------------------#
# Plot point estimates #
#----------------------#

# set position *pos* of each method
pos <- 0
repres <- cbind(repres, rep(0, dim(repres)[1]))
colnames(repres)[length(colnames(repres))] <- "pos"
for (meth in sort(unique(repres$method))) {
  repres[which(repres$method == meth), "pos"] <- pos
  pos <- pos + 1
}

# plot and save!!
point_estimates <- ggplot(data = repres,
                          aes(x = as.factor(s_true), y = est,
                              color = "replicates", fill  = "replicates")) +
  facet_grid(vars(method), vars(sample)) +
  geom_point(shape = 16, alpha = 0.7, size = 0.5,
             position = position_jitterdodge(jitter.width = 0.7)) +
  geom_boxplot(alpha = 0, col = "grey25", outlier.shape = NA, linewidth = 0.25) +
  geom_point(aes(x = as.factor(s_true), y = s_true),
             color = "#E41A1C", shape = 8, stroke = 0.25) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "grey") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "True s", y = "Estimated s", fill = "", color = "") +
  theme(plot.title = element_text(hjust = "0.5", face = "bold", size = 16),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(color = "grey95"),
        legend.text = element_text(size = 8),
        legend.direction = "vertical",
        legend.title = element_text(size = 10, vjust = .5, hjust = .3),
        legend.box = "vertical",
        legend.justification = "center",
        legend.box.just = "left",
        plot.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA),
        legend.position = "none",
        legend.key = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

point_estimates

ggsave(point_estimates,
       file = paste0(cwd, "../plots/panels/Figure_02/02-I-III_point_estimates.png"),
       width = 5.25, height = 7, dpi = 600, scale = 1.2)

#------------------------#
#          RMSE          #
#------------------------#

rmse <- ggplot(data = res[order(res$method), ],
               aes(x = s_true, y = rmse)) +
  facet_wrap(vars(sample_rev), nrow = 1, ncol = 3) +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = method, group = method),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = method),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = method, shape = method),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.1, 0.05),
                     minor_breaks = seq(0, 0.1, 0.01),
                     limits = c(0, 0.1)) +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt,
                     limits = c(0, 0.1)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1, 4, 3, 2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1, 4, 3, 2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 23, 21, 24),
                     name = "Method") +
  xlab("True s") +
  ylab("RMSE") +
  ggtitle("A") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.25)
  )

rmse

ggsave(rmse, file = paste0(cwd, "../plots/panels/Figure_03/03-I-III_a_rmse.pdf"),
       width = 5, height = 2, scale = 1.5)

#------------------------#
#          Bias          #
#------------------------#

bias <- ggplot(data = repres_stats[order(repres_stats$method), ],
               aes(x = s_true, y = Bias, col = as.factor(method))) +
  facet_wrap(vars(sample), nrow = 1, ncol = 3) +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt,
                     limits = c(0, 0.1)) +
  scale_y_continuous(breaks = seq(-150, 25, 50) / 1000,
                     minor_breaks = seq(-0.15, 0.025, 0.01),
                     limits = c(-0.1, 0.015)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1, 4, 3, 2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1, 4, 3, 2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 23, 21, 24),
                     name = "Method") +
  xlab("True s") +
  ylab("Bias") +
  ggtitle("B") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.25)
  )

bias

ggsave(bias, file = paste0(cwd, "../plots/panels/Figure_03/03-I-III_b_bias.pdf"),
       width = 5, height = 2, scale = 1.5)

#----------------------------#
#          Variance          #
#----------------------------#

variance <- ggplot(data = repres_stats[order(repres_stats$method), ],
                   aes(x = s_true, y = Variance, col = as.factor(method))) +
  facet_wrap(vars(sample), nrow = 1, ncol = 3) +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.001, 0.0005),
                     minor_breaks = seq(0, 0.001, 0.0001),
                     limits = c(0, 0.001)) +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1, 4, 3, 2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1, 4, 3, 2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 23, 21, 24),
                     name = "Method") +
  xlab("True s") +
  ylab("Variance") +
  ggtitle("C") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.25)
  )

variance

ggsave(variance, file = paste0(cwd, "../plots/panels/Figure_03/03-I-III_c_variance.pdf"),
       width = 5, height = 2, scale = 1.5)

#-------------------------#
# Plot True Positive Rate #
#-------------------------#

# since we cannot do this for bmws, make another dataframe without it
bmws <- which(res$method == "BMWS")
res_noBMWS <- res[-bmws, ]

# TPR doesn't make sense if true s is 0 (neutrality), so remove this case
neutral <- which(res_noBMWS$s_true == 0)
res_tpr <- res_noBMWS[-neutral, ]

# plot and save!!
tpr <- ggplot(data = res_tpr[order(res_tpr$method), ],
              aes(x = s_true, y = propPos, col = as.factor(method))) +
  facet_wrap(vars(sample_rev), nrow = 1, ncol = 3) +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x = -1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5),
                     minor_breaks = seq(0, 1.5, 0.1),
                     labels = paste0(seq(0, 150, 50), "%"),
                     limits = c(0, 1.02)) +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt,
                     limits = c(0, 0.1)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Accent")[c(1, 3, 2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 3, name = "Accent")[c(1, 3, 2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 21, 24),
                     name = "Method") +
  xlab("True s") +
  ylab("Power") +
  ggtitle("D") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.25)
  )

tpr

ggsave(tpr, file = paste0(cwd, "../plots/panels/Figure_03/03-I-III_d_tpr.pdf"),
       width = 5, height = 2, scale = 1.5)

# ApproxWF and Slattice are very close, plotting again and zooming onto
# y = 75% - 100% to see whether we notice a difference

tpr_zoom <- tpr + scale_y_continuous(breaks = seq(0, 1.5, 0.05),
                                     minor_breaks = seq(0, 1.5, 0.01),
                                     labels = paste0(seq(0, 150, 5), "%"),
                                     limits = c(0.75, 1.02))

tpr_zoom

ggsave(tpr_zoom,
       file = paste0(cwd, "../plots/panels/Figure_03/03-I-III_d_tpr_zoom.pdf"),
       width = 5, height = 2, scale = 1.5)

#--------------------------#
# Plot False Positive Rate #
#--------------------------#

# this time, we want to look only at the case where true s = 0

fpr <- ggplot(data = res_noBMWS[neutral, ],
              aes(x = as.factor(sample_rev), y = propPos)) +
  facet_wrap(vars(method), nrow = 1, ncol = 3) +
  geom_segment(x = -1, xend = 4, y = 0, yend = 0,
               color = "white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x = -1, xend = 4, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_bar(stat = "identity", fill = "grey", width = 0.65) +
  scale_y_continuous(breaks = seq(0, 0.09, 0.05),
                     minor_breaks = seq(0, 0.09, 0.01),
                     labels = paste0(seq(0, 9, 5), "%"),
                     limits = c(0, 0.06)) +
  xlab("Sampling Scheme") +
  ylab("False Positive Rate") +
  ggtitle("E") +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size = 8),
        legend.direction = "vertical",
        legend.title = element_text(size = 10, vjust = .5, hjust = .3),
        legend.box = "vertical",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        legend.key = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.y = element_line(color = "grey95"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.25)
  )

fpr

ggsave(fpr, file = paste0(cwd, "../plots/panels/Figure_03/03-I-III_e_fpr.pdf"), width = 2,
       height = 3, scale = 1.5)

#----------------------------------#
#       FINAL FIGURES 2 AND 3      #
#----------------------------------#

# removing legends
legend <- g_legend(rmse)
rmse <- rmse + theme(legend.position = "none")
bias <- bias + theme(legend.position = "none")
variance <- variance + theme(legend.position = "none")
tpr <- tpr + theme(legend.position = "none")

# removing facet grid strips on bias, variance, and tpr
bias <- bias + theme(strip.text = element_blank(),
                     strip.background = element_blank())
variance <- variance + theme(strip.text = element_blank(),
                             strip.background = element_blank())
tpr <- tpr + theme(strip.text = element_blank(),
                   strip.background = element_blank())

# removing x axis on rmse, bias, and variance
rmse <- rmse + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())
bias <- bias + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())
variance <- variance + theme(axis.title.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             axis.text.x = element_blank())

# first figure with only point-estimates
fig_I.III_1 <- point_estimates

fig_I.III_1

ggsave(fig_I.III_1, file = paste0(cwd, "../plots/Figure_02_ideal-dataset_point-estimates.png"),
       width = 6, height = 7, scale = 1.2)

# second figure with rmse, bias, variance, tpr, and fpr
fig_I.III_2 <- rmse / bias / variance / tpr / legend / fpr +
  plot_layout(heights = c(2.05, 2, 2,  2.05, 0.2, 1.7))

fig_I.III_2

ggsave(fig_I.III_2, file = paste0(cwd, "../plots/Figure_03_ideal-dataset-stats.pdf"),
       width = 6, height = 9, scale = 1.2)
