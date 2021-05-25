#!usr/bin/env R
# epitope_saliency.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads trained feature weights for models using epitope data and uses
# it to generate saliency plots for learned feature weights

# set proper working directory
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/plots")
# load in required libraries
library(ggplot2)
library(data.table)
library(reshape2)
library(ggseqlogo)
library(gridExtra)
library(cowplot)

# load in saliencies from motif model
motif <- fread("epitope_motif_saliency.csv", header = TRUE)
colnames(motif)[1:4] <- c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy")
motif$pos = rep(1:17, 2)

# subset to true cleavage examples
pos <- subset(motif[,c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy", "pos")], motif$cleave_type == "y")
norm_pos <- data.frame(apply(pos[,1:4], 2, scale, center=FALSE))
norm_pos$pos <- pos$pos

# reformat for plotting
pos_long <- melt(norm_pos, id.vars = "pos")
pos_long$variable <- gsub("[.]", " ", pos_long$variable)
pos_long$variable <- factor(pos_long$variable, levels = c("Conformational Entropy", "Hydrophobicity", "Molecular Volume", "Polarity"))
position_labels <- c(paste0("'", 8:1), as.character(0:8))

# plot
pos_motif <- ggplot(pos_long, aes(x=pos, y=variable, fill=value)) + 
  scale_fill_gradient(low = "black", high = "yellow")+
  geom_tile() + 
  scale_x_continuous(breaks=pos_long$pos, labels = rep(position_labels,4))+
  geom_vline(xintercept = 9.5, color = "white", lty=2) +
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


# load in sequence data
seq <- fread("epitope_sequence_saliency.csv", header = TRUE)
colnames(seq)[1:20] <- paste0("v", colnames(seq)[1:20])
# subset to cleavage exmaples
pos_seq <- subset(seq[,1:20], seq$cleave_type == "y")

# reformat to labeled matrix
seq_plot_data <- t(as.matrix(pos_seq[,1:20]))
aas <- c('A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
rownames(seq_plot_data) <- aas

# plot
pos_logo <- ggseqlogo(seq_plot_data) + 
  geom_vline(xintercept = 9.5, color = "black", lty=2) +
  theme(axis.text.x = element_blank(),
        # axis.text.y=element_blank(),
        # axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        # panel.grid.minor.y=element_blank(),
        # panel.grid.major.y=element_blank(),
        legend.position = "right",
        plot.margin = unit(c(1,1,0.2,1), "cm"))
pos_logo

# combine plots
out1 <- plot_grid(pos_logo, pos_motif, align = "v", rel_heights = c(2/3, 1/3), ncol=1, scale = 1.15, axis = "t")
out1

# output
ggsave2("~/PycharmProjects/pepsickle-paper/figures/supplement/epi_model_feature_saliency.png", out1, device = "png", width = 8, height = 5, units = "in")

