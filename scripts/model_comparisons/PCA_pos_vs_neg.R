#!usr/bin/env R

# PCA_pos_vs_neg.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# 

# setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/plots")
setwd("./pepsickle-paper/data/validation_data/plots")

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# dat <- fread("validation_features.csv", header = T)
dat <- fread("cleavage_training_windows.csv", header=T)
dat$id <- "digestion_training_data"

dat_26S <- fread("cleavage_26S_training_windows.csv", header = T)
dat_26S$id <- "digestion_26S_training_data"
dat <- rbind(dat, dat_26S)

ep_dat <- fread("epitope_training_windows.csv", header = T)
ep_dat$id <- "epitope_training_data"
dat <- rbind(dat, ep_dat)

bg <- fread("whole_proteome_frag_features.csv", header = T)
bg$class <- "bg"
dat <- rbind(dat, bg)

rm(list=c("bg", "ep_dat", "dat_26S"))

std_dat = prcomp(dat[,2:(dim(dat)[2]-2)], scale. = T)
loadings_1 <- std_dat$rotation[,1]
v_1 <- seq(1,52, 4)
v_2 <- seq(2,52, 4)
v_3 <- seq(3,52, 4)
v_4 <- seq(4,52, 4)

p1 <- qplot(x=seq(1,13), y=loadings_1[v_1])
p2 <- qplot(x=seq(1,13), y=loadings_1[v_2])
p3 <- qplot(x=seq(1,13), y=loadings_1[v_3])
p4 <- qplot(x=seq(1,13), y=loadings_1[v_4])
grid.arrange(p1,p2,p3,p4, ncol=2)


plot_vals = data.frame(std_dat$x[,1:10])
plot_vals$entry_id = 1:nrow(plot_vals)
plot_vals$source <- (dat$id)
plot_vals$group_id <- paste(dat$id, dat$class, sep="_")

rm(list=c("dat"))

colors = brewer.pal(length(unique(plot_vals$group_id)), "Paired")

plt <- ggplot(plot_vals[(plot_vals$source == 'digestion_26S_training_data'| 
                           plot_vals$source == 'digestion_training_data'),], 
              aes(x=PC1, y=PC2, colour=group_id)) + geom_point(alpha=.5)+
  scale_color_manual(values = colors)
plt

kplt2 <- ggplot(plot_vals[(plot_vals$source == 'digestion_26S_training_data'| 
                           plot_vals$source == 'digestion_training_data'),], 
              aes(x=PC3, y=PC2, colour=group_id)) + geom_point(alpha=.5)+
  scale_color_manual(values = colors)
plt2

grid.arrange(plt, plt2, nrow=1)


d1p <- ggplot(plot_vals[(plot_vals$group_id == 'digestion_training_data_pos'),], 
              aes(x=PC1, y=PC2))+geom_density_2d_filled()+guides(fill=FALSE)+ggtitle("20S pos")
d1n <- ggplot(plot_vals[(plot_vals$group_id == 'digestion_training_data_neg'),], 
              aes(x=PC1, y=PC2))+geom_density_2d_filled()+guides(fill=FALSE)+ggtitle("20S neg")

d2p <- ggplot(plot_vals[(plot_vals$group_id == 'digestion_26S_training_data_pos'),], 
              aes(x=PC1, y=PC2))+geom_density_2d_filled()+guides(fill=FALSE)+ggtitle("26S pos")
d2n <- ggplot(plot_vals[(plot_vals$group_id == 'digestion_26S_training_data_neg'),], 
              aes(x=PC1, y=PC2))+geom_density_2d_filled()+guides(fill=FALSE)+ggtitle("26S neg")

d3p <- ggplot(plot_vals[(plot_vals$group_id == 'epitope_training_data_pos'),], 
              aes(x=PC1, y=PC2))+geom_density_2d_filled()+guides(fill=FALSE)+ggtitle("epitope pos")
d3n <- ggplot(plot_vals[(plot_vals$group_id == 'epitope_training_data_neg'),], 
              aes(x=PC1, y=PC2))+geom_density_2d_filled()+guides(fill=FALSE)+ggtitle("epitope neg")

grid.arrange(d1p,d1n,d2p,d2n,d3p,d3n, ncol=2)

x_dens <- kde2d(plot_vals[(plot_vals$group_id == 'digestion_26S_training_data_pos'),1], 
                plot_vals[(plot_vals$group_id == 'digestion_26S_training_data_pos'),2])
y_dens <- kde2d(plot_vals[(plot_vals$group_id == 'digestion_training_data_pos'),1], 
                plot_vals[(plot_vals$group_id == 'digestion_training_data_pos'),2])
dense_diff <- x_dens$z - y_dens$z

plot_dat <- expand.grid(x_dens$x, x_dens$y)
plot_dat$density <- as.vector(t(dense_diff))

dense1 <- ggplot(plot_dat, aes(x=Var1, y=Var2, z=density))+
  stat_contour_filled()+
  xlim(min(plot_vals$PC1), max(plot_vals$PC1))+
  ylim(min(plot_vals$PC2), max(plot_vals$PC2))


x_dens2 <- kde2d(plot_vals[(plot_vals$group_id == 'digestion_26S_training_data_neg'),1], 
                plot_vals[(plot_vals$group_id == 'digestion_26S_training_data_neg'),2])
y_dens2 <- kde2d(plot_vals[(plot_vals$group_id == 'digestion_training_data_neg'),1], 
                plot_vals[(plot_vals$group_id == 'digestion_training_data_neg'),2])
dense_diff2 <- x_dens2$z - y_dens2$z

plot_dat2 <- expand.grid(x_dens2$x, x_dens2$y)
plot_dat2$density <- as.vector(t(dense_diff2))

dense2 <- ggplot(plot_dat2, aes(x=Var1, y=Var2, z=density))+
  stat_contour_filled()+
  xlim(min(plot_vals$PC1), max(plot_vals$PC1))+
  ylim(min(plot_vals$PC2), max(plot_vals$PC2))

grid.arrange(dense1, dense2, ncol=2)
