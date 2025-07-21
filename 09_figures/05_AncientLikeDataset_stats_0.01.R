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

repres <- read.table(paste0(cwd, "../results/04_II_ancient-like_allrepres.tsv"),
                     header = TRUE)

# do some cleaning up
repres[repres$method == "approxwf", "method"] <- "ApproxWF"
repres[repres$method == "bmws", "method"] <- "BMWS"
repres[repres$method == "slattice", "method"] <- "Slattice"
repres[repres$method == "sr", "method"] <- "Sr"

repres$method <- factor(repres$method, levels = c("ApproxWF",
                                                  "BMWS",
                                                  "Slattice",
                                                  "Sr"))

repres$init_old <- repres$init
repres$init <- paste0("IAF = ", repres$init)

#-------------------------------------------------#
#          Load summary stats into "res"          #
#-------------------------------------------------#

res <- read.table(paste0(cwd, "../results/04_II_ancient-like_allres.tsv"),
                  header = TRUE)

res[res$method == "approxwf", "method"] <- "ApproxWF"
res[res$method == "bmws", "method"] <- "BMWS"
res[res$method == "slattice", "method"] <- "Slattice"
res[res$method == "sr", "method"] <- "Sr"

# do some cleaning up
res$method <- factor(res$method, levels = c("ApproxWF",
                                            "BMWS",
                                            "Slattice",
                                            "Sr"))

# considering only initial allele frequencies of 1%
res <- res[which(res$init == 0.01), ]
repres <- repres[which(repres$init_old == 0.01), ]

#----------------------------------------------#
#          Computing additional stats          #
#----------------------------------------------#

repres_stats <- data.frame()
head(repres)
for (method in unique(repres$method)) {
  for (s_true in unique(repres$s_true)) {
    for (init in unique(repres$init)) {
      
      temp <- repres[which(repres$method == method), ]
      temp <- temp[which(temp$s_true == s_true), ]
      temp <- temp[which(temp$init == init), ]
      
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
                          init = init,
                          init_old = unique(temp$init_old),
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

head(repres_stats)

#------------------------#
#          RMSE          #
#------------------------#

rmse <- ggplot(data = res[order(res$method), ],
               aes(x=s_true, y=rmse, col = as.factor(method))) +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.1),
                     minor_breaks = seq(0, 0.8, 0.05),
                     limits = c(0, 0.8)) +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,4,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,4,3,2)],
                     name = "Method") +
  scale_shape_manual(values = c(22, 23, 21, 24),
                     name = "Method") +
  xlab("True s") +
  ylab("RMSE") +
  ggtitle("A") +
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
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

rmse

ggsave(rmse, file = paste0(cwd, "../plots/panels/Figure_05/05-IV-0.01_a_rmse.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#------------------------#
#          Bias          #
#------------------------#

bias <- ggplot(data = repres_stats[order(repres_stats$method), ],
               aes(x=s_true, y=bias, col = as.factor(method))) +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(-0.1, 0.15, 0.05),
                     minor_breaks = seq(-0.1, 0.15, 0.01)) +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,4,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,4,3,2)],
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
        legend.text = element_text(size=8),
        legend.direction = "horizontal",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

bias

ggsave(bias, file = paste0(cwd, "../plots/panels/Figure_05/05-IV-0.01_b_bias.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#----------------------------#
#          Variance          #
#----------------------------#

variance <- ggplot(data = repres_stats[order(repres_stats$method), ],
                   aes(x=s_true, y=var, col = as.factor(method))) +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1),
                     minor_breaks = seq(0, 0.6, 0.05),
                     limits = c(0, 0.6)) +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt) +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,4,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Accent")[c(1,4,3,2)],
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
        legend.text = element_text(size=8),
        legend.direction = "horizontal",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

variance

ggsave(variance, file = paste0(cwd, "../plots/panels/Figure_05/05-IV-0.01_c_variance.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#-------------------------#
# Plot True Positive Rate #
#-------------------------#

# removing bmws because we can't do TPR with it
bmws <- which(res$method == "BMWS")
res_noBMWS <- res[-bmws, ]

# removing neutral case
neutral <- which(res_noBMWS$s_true == 0)
res_noBMWS_tpr <- res_noBMWS[-neutral, ]

# plot and save!!
tpr <- ggplot(data = res_noBMWS_tpr[order(res_noBMWS_tpr$method), ],
              aes(x=s_true, y=propPos, col = as.factor(method))) +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 1, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_line(aes(col = as.factor(method), group = as.factor(method)),
            linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = as.factor(method)),
             size = 4, stroke = 0.5, fill = "white", color = "white") +
  geom_point(aes(fill = as.factor(method), shape = as.factor(method)),
             size = 3, stroke = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(0, 1, 0.5),
                     minor_breaks = seq(0, 1, 0.1),
                     labels = paste0(seq(0, 100, 50), "%")) +
  scale_x_continuous(breaks = c(0, seq(0.01, 0.1, 0.01)),
                     labels = s_vals_txt) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
                    name = "Method") +
  scale_color_manual(values = brewer.pal(n = 3, name = "Accent")[c(1,3,2)],
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
        legend.text = element_text(size=8),
        legend.direction = "horizontal",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", color=NA),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

tpr

ggsave(tpr, file = paste0(cwd, "../plots/panels/Figure_05/05-IV-0.01_d_tpr.pdf"),
       width = 2.5, height = 2, scale = 1.5)

#--------------------------#
# Plot False Positive Rate #
#--------------------------#

#keeping only neutral case
res_noBMWS_fpr <- res_noBMWS[neutral, ]

fpr <- ggplot(data = res_noBMWS_fpr[order(res_noBMWS_fpr$method), ],
              aes(x=as.factor(method), y=propPos)) +
  geom_segment(x=-1, xend = 4, y = 0, yend = 0,
               color="white", linewidth = 0.5, linetype = "solid") +
  geom_segment(x=-1, xend = 4, y = 0, yend = 0,
               color = "grey65", linewidth = 0.25, linetype = "solid") +
  geom_bar(stat="identity", fill = "grey", width = 0.65) +
  scale_y_continuous(breaks=seq(0, 0.1, 0.01),
                     minor_breaks = seq(0, 0.1, 0.005),
                     labels = paste0(seq(0, 10, 1), "%"),
                     limits = c(0, 0.06)) +
  xlab("Method") +
  ylab("False Positive Rate") +
  ggtitle("E") +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = -15, l = -5)),
        legend.text = element_text(size=8),
        legend.direction = "vertical",
        legend.title = element_text(size=10, vjust = .5, hjust = .3),
        legend.box = "vertical",
        legend.justification = "center",
        legend.box.just = "left",
        legend.box.background = element_rect(fill = "white", color=NA),
        legend.position = "right",
        legend.key = element_rect(fill = "white", color=NA),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        panel.grid.major.y = element_line(color = "grey95"),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color="black", linewidth = 0.25)
  )

fpr

ggsave(fpr, file = paste0(cwd, "../plots/panels/Figure_05/05-IV-0.01_e_fpr.pdf"), width = 2,
       height = 3, scale = 1.5)

#-----------------------------#
#       FINAL FIGURES 5       #
#-----------------------------#

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

# figure with rmse, bias, variance, tpr, and fpr
fig_IV_2 <- ((rmse / bias / variance + plot_layout(heights = c(1, 1, 1))) |
               (tpr / fpr + plot_layout(heights = c(1.075, 2)))) /
  legend + plot_layout(heights = c(9.8, 0.2))

fig_IV_2

ggsave(fig_IV_2, file = paste0(cwd, "../plots/Figure_05_ancient-dataset_0.01_stats.pdf"),
       width = 6, height = 6, scale = 1.2)
