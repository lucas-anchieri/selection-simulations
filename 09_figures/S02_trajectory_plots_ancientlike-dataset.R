library(data.table)
library(ggplot2)
library(RColorBrewer)

wd <- "../"
scratch <- "./data/simulation_Ne10000_Ngen1500"

AncientLike_sampling <- data.frame(gen = c(1, 27, 109, 111, 163, 181, 201, 221, 232,
                                           235, 296, 301, 303, 308, 311, 318, 319, 321,
                                           322, 341),
                                   n = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 3,
                                         4, 3, 4, 9, 50))

my_theme <- function() {
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
    plot.background = element_rect(fill='white', color=NA),
    legend.box.background = element_rect(fill = "white", color=NA),
    legend.key = element_rect(fill = "white", color=NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color="black", linewidth = 0.25),
    axis.title = element_text(size = 10)
  )
}

get_stats <- function(gen, data) {
  interval <- sort(unlist(data[gen, ]))[c(25, 975)]
  mean <- mean(unlist(data[gen, ]))
  stats <- data.frame(gen = gen,
                      mean = mean,
                      low = interval[1],
                      high = interval[2])
  return(stats)
}

all_traj <- data.frame()

all_init <- c(0.01, 0.1)

all_s <- c(0, 0.005, seq(0.01, 0.1, 0.01))

all_parameters <- data.frame(init = rep(all_init, each = 12),
                             s = rep(all_s, 2))

load_trajectories <- function(i, parameters) {
  s <- parameters[i, "s"]
  init <- parameters[i, "init"]

  traj <- fread(paste0(scratch, "init", as.character(init),
                       "/init", as.character(init), "_s", as.character(s),
                       "/init", as.character(init), "_s", as.character(s), ".tsv"))

  traj <- traj[1:341, 1:1000]

  traj <- traj / 20000

  stats <- as.data.frame(do.call(rbind,
                                 lapply(seq(341), get_stats, data = traj)))
  stats <- stats[order(stats$gen), ]
  stats$init <- init
  stats$s <- s

  stats$n <- 0
  stats[AncientLike_sampling$gen, "n"] <- AncientLike_sampling$n

  traj <- cbind(stats, traj)

  traj_melt <- as.data.frame(melt(as.data.table(traj),
                                  id.vars = c("gen", "init", "s",
                                              "mean", "low", "high", "n")))

  return(traj_melt)
}

all_traj <- do.call(rbind,
                    lapply(seq(24), load_trajectories, parameters = all_parameters))

all_traj$s_label <- factor(paste0("s = ", all_traj$s),
                           levels = paste0("s = ", c(seq(0.1, 0.01, -0.01), 0.005, 0)))

all_traj$init_label <- factor(paste0("IAF = ", all_traj$init),
                              levels = c("IAF = 0.01", "IAF = 0.1"))

head(all_traj)
dim(all_traj)
unique(all_traj$variable)

plot <- ggplot(data = all_traj) +
  facet_grid(s_label ~ init_label) +
  geom_line(aes(x = gen, y = value, group = variable),
            alpha = 0.1, color = "grey", linewidth = 0.1) +
  geom_point(data = all_traj[all_traj$gen %in% AncientLike_sampling$gen, ],
               aes(x = gen, y = mean),
               color = "black", shape = 4) +
  scale_y_continuous(breaks = c(0, 0.5, 1),
                     minor_breaks = c(seq(0, 1, 0.1)),
                     labels = c("0", "0.5", "1")) +
  xlab("Generations") +
  ylab("Allele Frequency") +
  my_theme()

# plot

ggsave(plot = plot, paste0(wd, "plots/Supplementary_Figure_S02_trajectories_ancient-dataset_trajectories.png"),
       width = 3, height = 9, unit = "in", dpi = 600,
       scale = 1.5)

