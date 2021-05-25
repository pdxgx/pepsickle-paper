#!usr/bin/env R

# gradients.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads and plots gradients of trained features from a trained DL model

# load required packages
setwd("~/PycharmProjects/proteasome/netchop/results")
library(ggplot2)
library(reshape2)
library(scales)
# library(gridExtra)
library(egg)

# load in gradients
dat <- read.csv("mean_grads.csv", stringsAsFactors = F)
dat_long <- melt(dat, id.vars = "pos")
dat_long$value <- abs(dat_long$value)

# plot for all values
p <- ggplot(dat_long, aes(x=pos, y=variable, fill=value)) + 
  scale_fill_gradient(low = "white", high = "red")+
  geom_tile() + 
  scale_x_continuous(breaks=dat_long$pos)+
  theme_classic()+
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "right")
p

# sum across positions
pos_importance <- apply(dat[,1:4], 1, function(x){sum(abs(x))})
property_importance <- apply(dat[,1:4], 2, function(x){sum(abs(x))})
property_names = factor(colnames(dat)[1:4], levels=colnames(dat)[1:4])

# plot position importance
p_pos <- ggplot(dat, aes(x=dat$pos, y = pos_importance))+ 
  geom_bar(stat = "identity", width = 1, colour="black", fill="lightgrey") + 
  theme_classic()+ 
  scale_y_continuous(position = "right")+
  theme(axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())
p_pos

# plot property importance
p_prop <- ggplot(as.data.frame(property_importance), aes(x=property_names, y = property_importance))+ 
  geom_bar(stat = "identity",width = 1, colour="black", fill="lightgrey") + theme_classic()+ coord_flip()+scale_y_reverse()
p_prop

# plot w/ void
p_void <- ggplot(as.data.frame(property_importance), aes(x=colnames(dat)[1:4], y = property_importance))+ geom_bar(stat = "identity", alpha=0) + theme_void()

# arrange
ggarrange(p_void, p_pos, p_prop, p, widths = c(1,2), padding = 0)

