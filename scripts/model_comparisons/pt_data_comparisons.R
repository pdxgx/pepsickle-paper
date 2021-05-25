
#!usr/bin/env R

# pt_data_comparisons.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads in pepsickle results from runs on patient data
# and provides summary/statistical comparisons

# set working directory
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/pt_data/pepsickle_pt_summaries")

# Load in required packages
library(ggplot2)
library(data.table)
library(gridExtra)
library(cowplot)

# load in data from Ott. et all
ott_dat <- fread("ott_pt_summary.csv")
ott_dat$response <- as.factor(ott_dat$response)
ott_dat$source <- "Ott et al."

# load in data from mupexi
mupexi_dat <- fread("mupexi_pt_summary.csv")
mupexi_dat$response <- as.factor(mupexi_dat$response)
mupexi_dat$source <- "MuPeXI"

# load in data from TESLA study
tesla_dat <- fread("tesla_pt_summary.csv")
tesla_dat$response <- as.factor(tesla_dat$response)
tesla_dat$source <- "TESLA"

# combine results from all three
full_dat <- rbind(ott_dat, tesla_dat, mupexi_dat)

# look at terminal cleavage by study
o_term <- ggplot(ott_dat, aes(x=term_prob, fill=response)) + geom_density(alpha=.4) + xlim(c(0,1))
o_term
m_term <- ggplot(mupexi_dat, aes(x=term_prob, fill=response)) + geom_density(alpha=.4)+ xlim(c(0,1))
m_term
t_term <- ggplot(tesla_dat, aes(x=term_prob, fill=response)) + geom_density(alpha=.4)+ xlim(c(0,1))
t_term

# plot terminal cleavage across all studies
full_term <- ggplot(full_dat, aes(x=term_prob, fill=response)) + 
  geom_density(alpha=.4)+ xlim(c(0,1))+
  xlab("C-terminal Cleavage Probability")
full_term

# plot terminal cleavage as boxplots
full_term_bp <- ggplot(full_dat, aes(x=response, y=term_prob, fill=response))+geom_boxplot(alpha=.4, fill="white")+
  geom_point(aes(color = source), position = position_jitter(width = .05), alpha=.6) +
  scale_x_discrete(labels = c("No", "Yes"))+
  scale_color_manual(values = c("green", "orange", "magenta"))+
  scale_fill_manual(values = c("darkgray", "darkgray"))+
  ylab("C-terminal Cleavage Probability")+
  xlab("Responsive")+
  theme_classic()+
  guides(fill=FALSE, color=FALSE)+
  theme(legend.text = element_text(size=12), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size=12))
full_term_bp

# facet by study ID
full_term_by_study <- full_term_bp+facet_wrap(~source)+ 
  theme(axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank())
full_term_by_study

# combine for export
out_plot <- plot_grid(full_term_bp, full_term_by_study, align = "h", rel_widths = c(1/3, 2/3), nrow=1, axis = "t")

# save
ggsave("~/PycharmProjects/pepsickle-paper/figures/main/pt_cleavage_predictions.png", out_plot, height = 5, width = 8, units = "in", device = "png")

## look at statistical associations
# Compare for full dataset
quantile(full_dat$term_prob) # look at upper quartile as cutoff
table(full_dat$term_prob >= 0.812475, full_dat$response) # 75th percentile
## PPV
# baseline
bl <- 45/(45+717)
# post-filter
pf <- 18/(18+173)
(pf-bl)/bl
18/(18+173)
45/(45+171)

chisq.test(table(full_dat$term_prob >= 0.812475, full_dat$response))

full_dat$term_cleaved2 <- 0
full_dat$term_cleaved2[full_dat$term_prob >= 0.81247] <- 1
caret::confusionMatrix(as.factor(full_dat$term_cleaved2), 
                       as.factor(full_dat$response), mode = "prec_recall", positive="1")


wilcox.test(log(full_dat$term_prob)~full_dat$response)
