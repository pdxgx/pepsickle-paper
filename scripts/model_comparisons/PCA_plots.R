#!usr/bin/env R

# PCA_plots.R
#
# For issues contact Ben Weeder (weeder@ohsu.edu)
#
# Generates PCA plots of the training window sampling space

# setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/plots")
# setwd("/home/groups/ThompsonLab/weederb/training_window_features")
setwd(".pepsickle-paper/data/validation_data/output/plots")

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(MASS)

# dat <- fread("validation_features.csv", header = T)
dat <- fread("20S_physical_training_windows21aa.csv", header=T)
dat$id <- "digestion_training_data"

# dat_26S <- fread("cleavage_26S_training_windows.csv", header = T)
# dat_26S$id <- "digestion_26S_training_data"
# dat <- rbind(dat, dat_26S)

ep_dat <- fread("epitope_physical_training_windows21aa.csv", header = T)
ep_dat$id <- "epitope_training_data"
dat <- rbind(dat, ep_dat, fill=T)

bg <- fread("whole_proteome_frag_features_21aa.csv", header = T)
bg$class <- "bg"
dat <- rbind(dat, bg, fill=T)

# rm(list=c("bg", "ep_dat", "dat_26S"))
rm(list=c("bg", "ep_dat"))

std_dat <- prcomp(dat[,2:(dim(dat)[2]-2)], scale. = T)

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

qplot(y=std_dat$sdev)
qplot(y=std_dat$rotation[,1])

rm(list=c("dat", "std_dat"))

colors <- brewer.pal(length(unique(plot_vals$group_id)), "Paired")
plt <- ggplot(plot_vals, aes(x=PC1, y=PC2, colour=group_id)) + geom_point()+scale_color_manual(values = colors)
plt

plt2 <- ggplot(plot_vals, aes(x=PC1, y=PC3, colour=group_id)) + geom_point()+scale_color_manual(values = colors)
plt2

plt3 <- ggplot(plot_vals, aes(x=PC3, y=PC2, colour=group_id)) + geom_point()+scale_color_manual(values = colors)
plt3

plt4 <- ggplot(plot_vals, aes(x=PC3, y=PC4, colour=group_id)) + geom_point()+scale_color_manual(values = colors)
plt4


d1 <- ggplot(plot_vals[plot_vals$coverage=="human_bg",], 
             aes(x=PC1, y=PC2))+
  geom_density2d_filled(n=200)+
  guides(fill=FALSE)+
  ylim(min(plot_vals$PC2), max(plot_vals$PC2))+
  xlim(min(plot_vals$PC1), max(plot_vals$PC1))+
  ggtitle("Human Background")

d2 <- ggplot(plot_vals[plot_vals$coverage=="epitope_sampled",], 
             aes(x=PC1, y=PC2))+
  geom_density2d_filled(n=200)+
  guides(fill=FALSE)+
  ylim(min(plot_vals$PC2), max(plot_vals$PC2))+
  xlim(min(plot_vals$PC1), max(plot_vals$PC1))+
  ggtitle("Epitope Training Data")

d3 <- ggplot(plot_vals[plot_vals$coverage=="cleavage_20S_sampled",], 
             aes(x=PC1, y=PC2))+
  geom_density2d_filled(n=200)+
  guides(fill=FALSE)+
  ylim(min(plot_vals$PC2), max(plot_vals$PC2))+
  xlim(min(plot_vals$PC1), max(plot_vals$PC1))+
  ggtitle("20S Training Data")

out_1 <- grid.arrange(d1,d2,d3, ncol=2)


bg_density <- kde2d(plot_vals[plot_vals$coverage=="human_bg", 1],
                    plot_vals[plot_vals$coverage=="human_bg", 2],
                    n=100, lims = c(min(plot_vals$PC1), max(plot_vals$PC1),
                                    min(plot_vals$PC2), max(plot_vals$PC2)))

bg_plot_dat <- expand.grid(bg_density$x, bg_density$y)
bg_plot_dat$density <- as.vector(t(bg_density$z))

c1 <- ggplot(bg_plot_dat, aes(x=Var2, y=Var1, z=density))+
  geom_contour_filled(breaks = c(0,1e-6,seq(1e-6,.1,.005)))+
  guides(fill=FALSE)+ggtitle("Human BG Data w/Footprint")

ep_density <- kde2d(plot_vals[plot_vals$coverage=="epitope_sampled",1],
                    plot_vals[plot_vals$coverage=="epitope_sampled",2],
                    n=100, lims = c(min(plot_vals$PC1), max(plot_vals$PC1),
                                    min(plot_vals$PC2), max(plot_vals$PC2)))

ep_plot_dat <- expand.grid(ep_density$x, ep_density$y)
ep_plot_dat$density <- as.vector(t(ep_density$z))

c2 <- ggplot(ep_plot_dat, aes(x=Var2, y=Var1, z=density))+
  geom_contour_filled(breaks = c(0,1e-6,seq(1e-6,.1,.005)))+
  guides(fill=FALSE)+ggtitle("Epitope Traning Data w/Footprint")

d20S_density <- kde2d(plot_vals[plot_vals$coverage=="cleavage_20S_sampled",1],
                    plot_vals[plot_vals$coverage=="cleavage_20S_sampled",2],
                    n=100, lims = c(min(plot_vals$PC1), max(plot_vals$PC1),
                                    min(plot_vals$PC2), max(plot_vals$PC2)))

d20S_plot_dat <- expand.grid(d20S_density$x, d20S_density$y)
d20S_plot_dat$density <- as.vector(t(d20S_density$z))

c3 <- ggplot(d20S_plot_dat, aes(x=Var2, y=Var1, z=density))+
  geom_contour_filled(breaks = c(0,1e-6,seq(1e-6,.1,.005)))+
  guides(fill=FALSE)+ggtitle("20S Traning Data w/Footprint")

out_2 <- grid.arrange(c1,c2, c3, ncol=1)


p0 <- ggplot(plot_vals, aes(x=PC1, y=PC2, colour=as.factor(coverage), group=as.factor(coverage)))+
  geom_point(data=plot_vals[plot_vals$coverage=="human_bg",], aes(x=PC1, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  scale_color_manual(values = c("grey")) + guides(colour=FALSE)

p1 <- ggplot(plot_vals, aes(x=PC1, y=PC2, colour=as.factor(coverage), group=as.factor(coverage)))+
  geom_point(data=plot_vals[plot_vals$coverage=="human_bg",], aes(x=PC1, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  geom_point(data=plot_vals[plot_vals$coverage=="epitope_sampled",], aes(x=PC1, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  scale_color_manual(values = c("purple", "grey")) + guides(colour=FALSE)

p2 <- ggplot(plot_vals, aes(x=PC3, y=PC2, colour=as.factor(coverage), group=as.factor(coverage)))+
  geom_point(data=plot_vals[plot_vals$coverage=="human_bg",], aes(x=PC3, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  geom_point(data=plot_vals[plot_vals$coverage=="epitope_sampled",], aes(x=PC3, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  scale_color_manual(values = c("purple", "grey"))+ guides(colour=FALSE)

p3 <- ggplot(plot_vals, aes(x=PC1, y=PC2, colour=as.factor(coverage), group=as.factor(coverage)))+
  geom_point(data=plot_vals[plot_vals$coverage=="human_bg",], aes(x=PC1, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  geom_point(data=plot_vals[plot_vals$coverage=="cleavage_20S_sampled",], aes(x=PC1, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  scale_color_manual(values = c("green", "grey")) + guides(colour=FALSE)

p4 <- ggplot(plot_vals, aes(x=PC3, y=PC2, colour=as.factor(coverage), group=as.factor(coverage)))+
  geom_point(data=plot_vals[plot_vals$coverage=="human_bg",], aes(x=PC3, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  geom_point(data=plot_vals[plot_vals$coverage=="cleavage_20S_sampled",], aes(x=PC3, y=PC2, colour=as.factor(coverage)), alpha=.5)+
  scale_color_manual(values = c("green", "grey"))+ guides(colour=FALSE)


p_all <- grid.arrange(p1,p2,p3,p4, ncol = 2)
print(p_all)

ggsave("out1.png", plot = out_1, device = "png")
ggsave("out2.png", plot = out_2, device = "png")
ggsave("out3.png", plot = p_all, device = "png")
