#!usr/bin/env R

# UMAP_plotting.R
#
# For issues contact Ben Weeder (weeder@ohsu.edu)
#
# This script generates UMAP projections of the sampled feature space

require(umap)
require(data.table)
require(ggplot2)

# setwd("/home/groups/ThompsonLab/weederb/training_window_features")
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/plots")

dat <- fread("20S_physical_training_windows21aa.csv", header=T)
dat$id <- "digestion_training_data"

ep_dat <- fread("epitope_physical_training_windows21aa.csv", header = T)
ep_dat$id <- "epitope_training_data"
dat <- rbind(dat, ep_dat)

bg <- fread("whole_proteome_frag_features_21aa.csv", header = T)
bg$class <- "bg"
dat <- rbind(dat, bg)

std_dat <- prcomp(dat[,2:(dim(dat)[2]-2)], scale. = T)

# generate plotable data
plot_vals <- data.frame(std_dat$x[,1:10])
plot_vals$entry_id <- 1:nrow(plot_vals)
plot_vals$group_id <- dat$id
plot_vals$coverage <- "cleavage_20S_sampled"
# plot_vals$coverage[plot_vals$group_id == "digestion_26S_training_data"] <- "cleavage_26S_sampled"
plot_vals$coverage[plot_vals$group_id == "epitope_training_data"] <- "epitope_sampled"
plot_vals$coverage[plot_vals$group_id == "human_proteome"] <- "human_bg"
# plot_vals$coverage <- factor(plot_vals$coverage, levels=c("human_bg", "cleavage_20S_sampled", "cleavage_26S_sampled", "epitope_sampled"))
plot_vals$coverage <- factor(plot_vals$coverage, levels=c("human_bg", "cleavage_20S_sampled", "epitope_sampled"))
levels(plot_vals$coverage)


umap_vals <- umap(plot_vals[,0:10])
print(str(umap_vals))
umap_p <- ggplot(data=umap_vals, aes(x=umap_vals[,1], y=umap_vals[,2]))+geom_point()
ggsave(filename = "out_umap_test.pdf", plot = umap_p, device = "png")
                 
