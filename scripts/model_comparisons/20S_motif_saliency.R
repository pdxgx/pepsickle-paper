#!usr/bin/env R

# 20S_sequence_saliency.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads trained feature weights for models using 20S data and uses
# it to generate saliency plots for learned feature weights

# Load in required packages
library(ggplot2)
library(data.table)
library(reshape2)
library(ggseqlogo)
library(gridExtra)

# set working directory to location of saliency vals
# (change wd to proper location)
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/plots")

# read in saliency values from motif model
motif <- fread("20S_motif_saliency.csv", header = TRUE)
# rename columns and add position info
colnames(motif)[1:4] <- c("v1", "v2", "v3", "v4")
motif$pos = rep(1:13, 4)

# subset to only saliencies for positive cleavage examples
c_pos <- subset(motif[,c("v1","v2","v3","v4","pos")], motif$cleave_type == "y" & motif$prot_type == "c")
c_norm_pos <- data.frame(apply(c_pos[,1:4], 2, scale, center=FALSE))
c_norm_pos$pos <- c_pos$pos

#reformat to long for plotting
c_pos_long <- melt(c_norm_pos, id.vars = "pos")

# generate saliency heatmap
p <- ggplot(c_pos_long, aes(x=pos, y=variable, fill=value)) + 
  scale_fill_gradient2(low = "blue", mid = "white", midpoint = 0, high = "red")+
  geom_tile() + 
  scale_x_continuous(breaks=c_pos_long$pos)+
  theme_classic()+
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "bottom")

p

# load in sequence saliency values
seq <- fread("20S_sequence_saliency.csv", header = TRUE)
# rename columns
colnames(seq)[1:20] <- paste0("v", colnames(seq)[1:20])
# subset to positive examples
c_pos_seq <- subset(seq[,1:20], seq$cleave_type == "y" & seq$prot_type == "c")

# reformat into named matrix with AA values
c_seq_plot_data <- t(as.matrix(c_pos_seq[,1:20]))
aas <- c('A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
rownames(c_seq_plot_data) <- aas

# generate logo plot from Sequence saliencies
p2 <- ggseqlogo(c_seq_plot_data, method='custom', seq_type="aa") + 
  geom_hline(size=1.5, yintercept=0, linetype=1)+
  ylim(-.2,.2)+
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "none")
p2

# generate plot with both sub-plots
grid.arrange(p2,p, ncol=1)

