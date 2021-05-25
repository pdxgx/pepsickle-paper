#!usr/bin/env R
# logo_plots.R
#
# For issues contact Ben Weeder (weeder@ohsu.edu)
#
# This script loads training windows for in-vitro and epitope data for
# subsequent visualization

#### logo plots comparing epitope and 20S training data
library(ggplot2)
library(data.table)
library(reshape2)
library(ggseqlogo)
library(gridExtra)
library(cowplot)

# setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/plots")
setwd("./pepsickle-paper/data/validation_data/output/plots")

cleavage_dat <- fread("20S_sequence_training_windows21aa.csv", header = T)
cleavage_dat <- cleavage_dat[,2:dim(cleavage_dat)[2]]
cleavage_chem_dat <- fread("20S_physical_training_windows21aa.csv", header=T)

epitope_dat <- fread("epitope_sequence_training_windows21aa.csv", header = T)
epitope_dat <- epitope_dat[,2:dim(epitope_dat)[2]]
epitope_chem_dat <- fread("epitope_physical_training_windows21aa.csv", header = T)


c_pos <- subset(cleavage_dat, cleavage_dat$class=="pos")
c_pos_freqs <- apply(c_pos[,1:(dim(c_pos)[2]-1)], 2, sum)
positions <- length(c_pos_freqs)/20
c_pos_formatted <- matrix(c_pos_freqs, nrow = 20, ncol = positions, byrow = T)
c_pos_formatted <- t(apply(c_pos_formatted, 1, function(x){x/apply(c_pos_formatted, 2, sum)}))
# c_pos_expanded <- cbind(matrix(data=0, nrow=nrow(c_pos_formatted), ncol=5), 
                        # c_pos_formatted,
                        # matrix(data=0, nrow=nrow(c_pos_formatted), ncol=5))

aas <- c('A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
rownames(c_pos_formatted) <- aas
# rownames(c_pos_expanded) <- aas

c_pos_logo <- ggseqlogo(c_pos_formatted, seq_type="aa") + 
  theme(axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,0,1), "cm"))
c_pos_logo


c_neg <- subset(cleavage_dat, cleavage_dat$class=="neg")
c_neg_freqs <- apply(c_neg[,1:(dim(c_neg)[2]-1)], 2, sum)
positions <- length(c_neg_freqs)/20
c_neg_formatted <- matrix(c_neg_freqs, nrow = 20, ncol = positions, byrow = T)
c_neg_formatted <- t(apply(c_neg_formatted, 1, function(x){x/apply(c_neg_formatted, 2, sum)}))
# c_neg_expanded <- cbind(matrix(data=0, nrow=nrow(c_neg_formatted), ncol=5), 
                        # c_neg_formatted,
                        # matrix(data=0, nrow=nrow(c_neg_formatted), ncol=5))

rownames(c_neg_formatted) <- aas
# rownames(c_neg_expanded) <- aas

c_neg_logo <- ggseqlogo(c_neg_formatted, seq_type="aa") + 
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "none")
c_neg_logo


grid.arrange(c_pos_logo, c_neg_logo)


e_pos <- subset(epitope_dat, epitope_dat$class=="pos")
e_pos_freqs <- apply(e_pos[,1:(dim(e_pos)[2]-1)], 2, sum)
positions <- length(e_pos_freqs)/20
e_pos_formatted <- matrix(e_pos_freqs, nrow = 20, ncol = positions, byrow = T)
e_pos_formatted <- t(apply(e_pos_formatted, 1, function(x){x/apply(e_pos_formatted, 2, sum)}))
# e_pos_expanded <- cbind(matrix(data=0, nrow=nrow(e_pos_formatted), ncol=5), 
# e_pos_formatted,
# matrix(data=0, nrow=nrow(e_pos_formatted), ncol=5))

aas <- c('A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
rownames(e_pos_formatted) <- aas
# rownames(e_pos_expanded) <- aas

e_pos_logo <- ggseqlogo(e_pos_formatted, seq_type="aa") + 
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,0,1), "cm"))
e_pos_logo


e_neg <- subset(epitope_dat, epitope_dat$class=="neg")
e_neg_freqs <- apply(e_neg[,1:(dim(e_neg)[2]-1)], 2, sum)
positions <- length(e_neg_freqs)/20
e_neg_formatted <- matrix(e_neg_freqs, nrow = 20, ncol = positions, byrow = T)
e_neg_formatted <- t(apply(e_neg_formatted, 1, function(x){x/apply(e_neg_formatted, 2, sum)}))
# e_neg_expanded <- cbind(matrix(data=0, nrow=nrow(e_neg_formatted), ncol=5), 
# e_neg_formatted,
# matrix(data=0, nrow=nrow(e_neg_formatted), ncol=5))

rownames(e_neg_formatted) <- aas
# rownames(e_neg_expanded) <- aas

e_neg_logo <- ggseqlogo(e_neg_formatted, seq_type="aa") + 
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "none")
e_neg_logo


grid.arrange(e_pos_logo, e_neg_logo)


grid.arrange(c_pos_logo, c_neg_logo, e_pos_logo, e_neg_logo, ncol=2)




cleavage_chem_dat <- fread("20S_physical_training_windows21aa.csv", header=T)
cleavage_chem_pos <- cleavage_chem_dat[cleavage_chem_dat$class == "pos",]
cleavage_chem_pos_values <- apply(cleavage_chem_pos[,2:(dim(cleavage_chem_pos)[2]-1)], 2, function(x){mean(abs(x))})
chem_pos_mat <- matrix(cleavage_chem_pos_values, ncol = 4, byrow=T)
chem_pos_mat <- as.data.frame(scale(chem_pos_mat))
colnames(chem_pos_mat) <- c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy")
chem_pos_mat$pos <- 1:nrow(chem_pos_mat)
chem_pos_long <- melt(chem_pos_mat, id.vars = "pos")
pos_labels <- c(paste0("'", 10:1), as.character(0:10))

chem_value_plot <- ggplot(chem_pos_long, aes(x=pos, y=variable, fill=value)) + 
  scale_fill_gradient(low = "black", high = "yellow")+
  geom_tile() + 
  scale_x_continuous(breaks=chem_pos_long$pos, labels = rep(pos_labels,4))+
  theme_classic()+
  xlab("Distance From Cleavage Point")+
  theme(axis.line.y=element_blank(),
        # axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,3,1,3), "cm"))
chem_value_plot

plot_grid(c_pos_logo, chem_value_plot, align = "v", rel_heights = c(2/3, 1/3), ncol=1, scale = 1.15)

# epitope
epi_chem_dat <- fread("epitope_physical_training_windows21aa.csv", header=T)
epi_chem_pos <- epi_chem_dat[epi_chem_dat$class == "pos",]
epi_chem_pos_values <- apply(epi_chem_pos[,2:(dim(epi_chem_pos)[2]-1)], 2, function(x){mean(abs(x))})
epi_chem_pos_mat <- matrix(epi_chem_pos_values, ncol = 4, byrow=T)
epi_chem_pos_mat <- as.data.frame(scale(epi_chem_pos_mat))
colnames(epi_chem_pos_mat) <- c("Polarity", "Molecular Volume", "Hydrophobicity", "Conformational Entropy")
epi_chem_pos_mat$pos <- 1:nrow(epi_chem_pos_mat)
epi_chem_pos_long <- melt(epi_chem_pos_mat, id.vars = "pos")
pos_labels <- c(paste0("'", 10:1), as.character(0:10))

epi_chem_value_plot <- ggplot(epi_chem_pos_long, aes(x=pos, y=variable, fill=value)) + 
  scale_fill_gradient(low = "black", high = "yellow")+
  geom_tile() + 
  scale_x_continuous(breaks=epi_chem_pos_long$pos, labels = rep(pos_labels,4))+
  theme_classic()+
  xlab("Distance From Cleavage Point")+
  theme(axis.line.y=element_blank(),
        # axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,3,1,3), "cm"))
epi_chem_value_plot

plot_grid(e_pos_logo, epi_chem_value_plot, align = "v", rel_heights = c(2/3, 1/3), ncol=1, scale = 1.15)

