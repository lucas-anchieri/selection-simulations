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

head(repres)

# plot and save!!
point_estimates <- ggplot(data = repres,
                          aes(x=as.factor(s_true), y=est,
                              color="replicates", fill ="replicates")) +
  facet_grid(vars(init), vars(method)) +
  geom_point(shape = 16, alpha=0.7, size= 0.5, position = position_jitterdodge(jitter.width = 0.7)) +
  geom_boxplot(alpha = 0, col = "grey25", outlier.shape = NA, linewidth = 0.25) +
  geom_point(aes(x=as.factor(s_true), y = s_true), color="#E41A1C", shape = 8, stroke = 0.25) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "grey") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "True s", y = "Estimated s", fill = "", color = "") +
  ylim(-0.025, 0.2) + # applying a limit because Sr goes way too high
  theme(plot.title = element_text(hjust = "0.5", face = "bold", size = 16),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_line(color = "grey95"),
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
        axis.line = element_line(color="black", linewidth = 0.25),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
  )

point_estimates

# save as separate panel
ggsave(point_estimates,
       file = paste0(cwd, "../plots/panels/Figure_04/04-IV_point_estimates.png"),
       width = 6, height = 3.7, dpi = 600, scale = 1.2)

#------------------------------#
#       FINAL FIGURES 4        #
#------------------------------#

# save as final figure
ggsave(point_estimates, file = paste0(cwd, "../plots/Figure_04_ancient-dataset_point-estimates.png"),
       width = 6, height = 3.7, scale = 1.2, dpi = 600)

