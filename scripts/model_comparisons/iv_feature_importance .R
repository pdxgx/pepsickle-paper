#!usr/bin/env R
# iv_feature_importance.R
#
# For issues contact Ben Weeder (weeder@ohsu.edu)
#
# This script plots the learned feature importances from the trained in-vitro ML
# model for visualization

# setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output")
setwd("./pepsickle-paper/data/validation_data/output/plots")

library(ggplot2)
library(data.table)
library(reshape2)
library(ggseqlogo)
library(gridExtra)
library(cowplot)

motif <- fread("final_feature_importances.tsv", header = TRUE)[1:28,1]
motif$variable <- rep(c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy"), 7)
motif$position <- rep(1:7, each=4)
position_labels <- c(paste0("'", 3:1), as.character(0:3))

pos_motif <- ggplot(motif, aes(x=position, y=variable, fill=mean_importance)) + 
  scale_fill_gradient(low = "black", high = "yellow")+
  geom_tile() + 
  geom_vline(xintercept = 4.5, color = "white", lty=2) +
  scale_x_continuous(breaks=motif$pos, labels = rep(position_labels, each=4))+
  theme_classic()+
  xlab("Distance From Cleavage Site")+
  labs(fill="Importance")+
  theme(axis.line.y=element_blank(), 
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size=12),
        # axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "right",
        plot.margin = unit(c(0,3,1,3), "cm"))
pos_motif

ggsave2("./plots/iv_model_feature_importance.png", pos_motif, device = "png", width = 8, height = 2, units = "in")




