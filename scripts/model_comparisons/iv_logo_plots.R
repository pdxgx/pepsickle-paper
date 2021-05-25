#!usr/bin/env R
# iv_logo_plots.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads in raw training window features for 20S fragments and uses
# it to generate feature frequency plots

# load in required libraries
library(ggplot2)
library(data.table)
library(reshape2)
library(ggseqlogo)
library(gridExtra)
library(cowplot)

## load in cleavage windows data
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/plots")
cleavage_dat <- fread("20S_sequence_training_windows7aa.csv", header = T)
cleavage_dat <- cleavage_dat[,2:dim(cleavage_dat)[2]]
cleavage_chem_dat <- fread("20S_physical_training_windows7aa.csv", header=T)

## create logo plots
# subsed to positive examples
c_pos <- subset(cleavage_dat, cleavage_dat$class=="pos")
c_pos_freqs <- apply(c_pos[,1:(dim(c_pos)[2]-1)], 2, sum)
positions <- length(c_pos_freqs)/20
c_pos_formatted <- matrix(c_pos_freqs, nrow = 20, ncol = positions, byrow = T)
c_pos_formatted <- t(apply(c_pos_formatted, 1, function(x){x/apply(c_pos_formatted, 2, sum)}))


aas <- c('A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
rownames(c_pos_formatted) <- aas

c_pos_logo <- ggseqlogo(c_pos_formatted, seq_type="aa") + 
  ylim(0,.2)+
  geom_vline(xintercept = 4.5, lty=2)+
  theme(axis.text.x = element_blank(),
        # axis.text.y=element_blank(),
        # axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        # panel.grid.minor.y=element_blank(),
        # panel.grid.major.y=element_blank(),
        legend.position = "right",
        plot.margin = unit(c(1,1,0.1,1), "cm"))
c_pos_logo

# repeat for negatives
c_neg <- subset(cleavage_dat, cleavage_dat$class=="neg")
c_neg_freqs <- apply(c_neg[,1:(dim(c_neg)[2]-1)], 2, sum)
positions <- length(c_neg_freqs)/20
c_neg_formatted <- matrix(c_neg_freqs, nrow = 20, ncol = positions, byrow = T)
c_neg_formatted <- t(apply(c_neg_formatted, 1, function(x){x/apply(c_neg_formatted, 2, sum)}))

rownames(c_neg_formatted) <- aas

c_neg_logo <- ggseqlogo(c_neg_formatted, seq_type="aa") + 
  geom_vline(xintercept = 4.5, lty=2)+
  ylim(0,.2)+
  theme(axis.text.x = element_blank(),
        # axis.text.y=element_blank(),
        # axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        # panel.grid.minor.y=element_blank(),
        # panel.grid.major.y=element_blank(),
        legend.position = "right",
        plot.margin = unit(c(1,1,0.1,1), "cm"))
c_neg_logo

## create chem heatmaps
# positive examples
cleavage_chem_pos <- cleavage_chem_dat[cleavage_chem_dat$class == "pos",]
cleavage_chem_pos_values <- apply(cleavage_chem_pos[,2:(dim(cleavage_chem_pos)[2]-1)], 2, function(x){mean(abs(x))})
chem_pos_mat <- matrix(cleavage_chem_pos_values, ncol = 4, byrow=T)
chem_pos_mat <- as.data.frame(scale(chem_pos_mat))
colnames(chem_pos_mat) <- c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy")
chem_pos_mat$pos <- 1:nrow(chem_pos_mat)
chem_pos_long <- melt(chem_pos_mat, id.vars = "pos")
pos_labels <- c(paste0("'", 3:1), as.character(0:3))

chem_value_plot <- ggplot(chem_pos_long, aes(x=pos, y=variable, fill=value)) + 
  scale_fill_gradient(low = "black", high = "yellow", limits=c(-2.1, 2.1))+
  geom_tile() + 
  geom_vline(xintercept = 4.5, lty=2, color="white")+
  scale_x_continuous(breaks=chem_pos_long$pos, labels = rep(pos_labels,4))+
  theme_classic()+
  xlab("Distance From Cleavage Site")+
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
chem_value_plot

out_1 <- plot_grid(c_pos_logo, chem_value_plot, align = "v", rel_heights = c(2/3, 1/3), ncol=1, scale = 1.15)

# negative examples
cleavage_chem_neg <- cleavage_chem_dat[cleavage_chem_dat$class == "neg",]
cleavage_chem_neg_values <- apply(cleavage_chem_neg[,2:(dim(cleavage_chem_neg)[2]-1)], 2, function(x){mean(abs(x))})
chem_neg_mat <- matrix(cleavage_chem_neg_values, ncol = 4, byrow=T)
chem_neg_mat <- as.data.frame(scale(chem_neg_mat))
colnames(chem_neg_mat) <- c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy")
chem_neg_mat$neg <- 1:nrow(chem_neg_mat)
chem_neg_long <- melt(chem_neg_mat, id.vars = "neg")
neg_labels <- c(paste0("'", 3:1), as.character(0:3))

chem_neg_value_plot <- ggplot(chem_neg_long, aes(x=neg, y=variable, fill=value)) + 
  scale_fill_gradient(low = "black", high = "yellow", limits=c(-2.1, 2.1))+
  geom_tile() + 
  geom_vline(xintercept = 4.5, lty=2, color="white")+
  scale_x_continuous(breaks=chem_neg_long$neg, labels = rep(neg_labels,4))+
  theme_classic()+
  xlab("Distance From Cleavage Site")+
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
chem_neg_value_plot

out_2 <- plot_grid(c_neg_logo, chem_neg_value_plot, align = "v", rel_heights = c(2/3, 1/3), ncol=1, scale = 1.15)

out_3 <- plot_grid(out_1, out_2, align = "v", ncol = 1)
out_3

# save plots
ggsave2("iv_positive_training_features.png", out_1, device = "png")
ggsave2("iv_negatvie_training_features.png", out_2, device = "png")
ggsave2("~/PycharmProjects/pepsickle-paper/figures/supplement/iv_full_training_features.png", out_3, device = "png", width = 8, height = 10)

