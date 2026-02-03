library(data.table)
library(ggplot2)
library(RColorBrewer)

wd <- "../"
scratch <- "./data/simulation_Ne10000_Ngen1500"

sampling_100gens <- data.frame(gen = seq(10, 100, 10),
                               n = rep(100, 10),
                               dataset = rep("100 Gen.", 10))
sampling_500gens <- data.frame(gen = seq(50, 500, 50),
                               n = rep(100, 10),
                               dataset = rep("500 Gen.", 10))
sampling_1000gens <- data.frame(gen = seq(100, 1000, 100),
                                n = rep(100, 10),
                                dataset = rep("1000 Gen.", 10))

Ideal_sampling <- rbind(sampling_100gens,
                        sampling_500gens,
                        sampling_1000gens)

Ideal_sampling$dataset <- factor(Ideal_sampling$dataset,
                                 levels = c("100 Gen.",
                                            "500 Gen.",
                                            "1000 Gen."))

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

all_s <- c(0, 0.005, seq(0.01, 0.1, 0.01))

load_trajectories <- function(i, parameters) {
  s <- parameters[i]
  init <- 0.1

  traj <- fread(paste0(scratch, "init", as.character(init),
                       "/init", as.character(init), "_s", as.character(s),
                       "/init", as.character(init), "_s", as.character(s), ".tsv"))

  traj <- traj[1:1000, 1:1000]

  traj <- traj / 20000

  stats <- as.data.frame(do.call(rbind,
                                 lapply(seq(1000), get_stats, data = traj)))
  stats <- stats[order(stats$gen), ]
  stats$init <- init
  stats$s <- s

  traj <- cbind(stats, traj)

  traj_melt <- as.data.frame(melt(as.data.table(traj),
                                  id.vars = c("gen", "init", "s",
                                              "mean", "low", "high")))

  return(traj_melt)
}

all_traj_1000 <- do.call(rbind,
                         lapply(seq(12), load_trajectories, parameters = all_s))

all_traj_1000$s_label <- factor(paste0("s = ", all_traj_1000$s),
                                levels = paste0("s = ", c(seq(0.1, 0.01, -0.01), 0.005, 0)))

all_traj_100 <- all_traj_1000[all_traj_1000$gen <= 100, ]
all_traj_100$variable <- paste0("100_", all_traj_100$variable)
all_traj_100$dataset <- "100 Gen."
all_traj_500 <- all_traj_1000[all_traj_1000$gen <= 500, ]
all_traj_500$variable <- paste0("500_", all_traj_500$variable)
all_traj_500$dataset <- "500 Gen."
all_traj_1000$variable <- paste0("1000_", all_traj_1000$variable)
all_traj_1000$dataset <- "1000 Gen."

all_traj <- rbind(all_traj_100, all_traj_500, all_traj_1000)

all_traj$dataset <- factor(all_traj$dataset,
                           levels = c("100 Gen.", "500 Gen.", "1000 Gen."))

head(all_traj)
dim(all_traj)
unique(all_traj$variable)

sampling_100 <- which(all_traj$gen %in% Ideal_sampling[Ideal_sampling$dataset == "100 Gen.", "gen"] & all_traj$dataset == "100 Gen.")
sampling_500 <- which(all_traj$gen %in% Ideal_sampling[Ideal_sampling$dataset == "500 Gen.", "gen"] & all_traj$dataset == "500 Gen.")
sampling_1000 <- which(all_traj$gen %in% Ideal_sampling[Ideal_sampling$dataset == "1000 Gen.", "gen"] & all_traj$dataset == "1000 Gen.")

sampling_gens <- c(sampling_100, sampling_500, sampling_1000)

plot <- ggplot(data = all_traj) +
  facet_grid(s_label ~ dataset, scales = "free_x") +
  geom_line(aes(x = gen, y = value, group = variable),
            alpha = 0.1, color = "grey", linewidth = 0.1) +
  geom_point(data = all_traj[sampling_gens, ],
               aes(x = gen, y = mean),
               color = "black", shape = 4) +
  scale_y_continuous(breaks = c(0, 0.5, 1),
                     minor_breaks = c(seq(0, 1, 0.1)),
                     labels = c("0", "0.5", "1")) +
  xlab("Generations") +
  ylab("Allele Frequency") +
  my_theme()

# plot

ggsave(plot = plot, paste0(wd, "plots/Supplementary_Figure_S01_ideal-dataset_trajectories.png"),
       width = 4.5, height = 9, unit = "in", dpi = 600,
       scale = 1.5)

