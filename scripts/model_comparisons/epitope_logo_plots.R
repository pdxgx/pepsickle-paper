#!usr/bin/env R

# epitope_logo_plots.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads epitope training window features and 
# generates plots showing frequencies

# load in required data
library(ggplot2)
library(data.table)
library(reshape2)
library(ggseqlogo)
library(gridExtra)
library(cowplot)

## load in epitope sequence data
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/plots")
epitope_dat <- fread("epitope_sequence_training_windows17aa.csv", header = T)
epitope_dat <- epitope_dat[,2:dim(epitope_dat)[2]]

# subset to positive data
e_pos <- subset(epitope_dat, epitope_dat$class=="pos")
e_pos_freqs <- apply(e_pos[,1:(dim(e_pos)[2]-1)], 2, sum)
positions <- length(e_pos_freqs)/20
e_pos_formatted <- matrix(e_pos_freqs, nrow = 20, ncol = positions, byrow = T)
e_pos_formatted <- t(apply(e_pos_formatted, 1, function(x){x/apply(e_pos_formatted, 2, sum)}))
aas <- c('A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
rownames(e_pos_formatted) <- aas

# plot
e_pos_logo <- ggseqlogo(e_pos_formatted, seq_type="aa") + 
  geom_vline(xintercept = 9.5, lty=2)+
  ylim(0,.5)+
  theme(axis.text.x = element_blank(),
        # axis.text.y=element_blank(),
        # axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        # panel.grid.minor.y=element_blank(),
        # panel.grid.major.y=element_blank(),
        legend.position = "right",
        plot.margin = unit(c(1,1,0.1,1), "cm"))
e_pos_logo

rm(list = c("e_pos", "e_pos_formatted"))
e_neg <- subset(epitope_dat, epitope_dat$class=="neg")
rm(list = c("epitope_dat"))

# repeat with non-cleavage windows
e_neg_freqs <- apply(e_neg[,1:(dim(e_neg)[2]-1)], 2, sum)
positions <- length(e_neg_freqs)/20
e_neg_formatted <- matrix(e_neg_freqs, nrow = 20, ncol = positions, byrow = T)
e_neg_formatted <- t(apply(e_neg_formatted, 1, function(x){x/apply(e_neg_formatted, 2, sum)}))
rownames(e_neg_formatted) <- aas

# plot
e_neg_logo <- ggseqlogo(e_neg_formatted, seq_type="aa") + 
  geom_vline(xintercept = 9.5, lty=2)+
  ylim(0,.5)+
  theme(axis.text.x = element_blank(),
        # axis.text.y=element_blank(),
        # axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        # panel.grid.minor.y=element_blank(),
        # panel.grid.major.y=element_blank(),
        legend.position = "right",
        plot.margin = unit(c(1,1,0.1,1), "cm"))
e_neg_logo

rm(list = c("e_neg", "e_neg_formatted"))

# repeat with chemical features (positives)
epi_chem_dat <- fread("epitope_physical_training_windows17aa.csv", header=T)
epi_chem_pos <- epi_chem_dat[epi_chem_dat$class == "pos",]
epi_chem_pos_values <- apply(epi_chem_pos[,2:(dim(epi_chem_pos)[2]-1)], 2, function(x){mean(abs(x))})
epi_chem_pos_mat <- matrix(epi_chem_pos_values, ncol = 4, byrow=T)
epi_chem_pos_mat <- as.data.frame(scale(epi_chem_pos_mat))
colnames(epi_chem_pos_mat) <- c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy")
epi_chem_pos_mat$pos <- 1:nrow(epi_chem_pos_mat)
epi_chem_pos_long <- melt(epi_chem_pos_mat, id.vars = "pos")
pos_labels <- c(paste0("'", 8:1), as.character(0:8))

# plot
epi_chem_value_plot <- ggplot(epi_chem_pos_long, aes(x=pos, y=variable, fill=value)) + 
  scale_fill_gradient(low = "black", high = "yellow", limits=c(-3,4))+
  geom_tile() + 
  geom_vline(xintercept = 9.5, lty=2, color="white")+
  scale_x_continuous(breaks=epi_chem_pos_long$pos, labels = rep(pos_labels,4))+
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
epi_chem_value_plot

out_1 <- plot_grid(e_pos_logo, epi_chem_value_plot, align = "v", rel_heights = c(2/3, 1/3), ncol=1, scale = 1.15)
out_1

rm(list=c("epi_chem_pos", "epi_chem_pos_long", "epi_chem_pos_mat"))

# repeat with chemical properties (negatives)
epi_chem_neg <- epi_chem_dat[epi_chem_dat$class == "neg",]
epi_chem_neg_values <- apply(epi_chem_neg[,2:(dim(epi_chem_neg)[2]-1)], 2, function(x){mean(abs(x))})
epi_chem_neg_mat <- matrix(epi_chem_neg_values, ncol = 4, byrow=T)
epi_chem_neg_mat <- as.data.frame(scale(epi_chem_neg_mat))
colnames(epi_chem_neg_mat) <- c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy")
epi_chem_neg_mat$neg <- 1:nrow(epi_chem_neg_mat)
epi_chem_neg_long <- melt(epi_chem_neg_mat, id.vars = "neg")
neg_labels <- c(paste0("'", 8:1), as.character(0:8))

# plot
epi_neg_chem_value_plot <- ggplot(epi_chem_neg_long, aes(x=neg, y=variable, fill=value)) + 
  scale_fill_gradient(low = "black", high = "yellow", limits=c(-3,4))+
  geom_tile() + 
  geom_vline(xintercept = 9.5, lty=2, color="white")+
  scale_x_continuous(breaks=epi_chem_neg_long$neg, labels = rep(neg_labels,4))+
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
# epi_neg_chem_value_plot

out_2 <- plot_grid(e_neg_logo, epi_neg_chem_value_plot, align = "v", rel_heights = c(2/3, 1/3), ncol=1, scale = 1.15)

out_1
out_2

out_3 <- plot_grid(out_1, out_2, align = "v", ncol = 1)

# save plots
ggsave2("epi_positive_training_features.png", out_1, device = "png")
ggsave2("eip_negatvie_training_features.png", out_2, device = "png")
ggsave2("~/PycharmProjects/pepsickle-paper/figures/supplement/epitope_full_training_features.png", out_3, device = "png", width = 8, height = 10)
