
#!usr/bin/env R

# epitope_mod_window_comparison.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads predictions from varying test window sizes and compares performance.

# Load in required packages
# setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/epitope_mod_window_preds")
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/epitope_mod_window_preds")

library(pROC)
library(ggplot2)
library(data.table)

# load DL predictions from each window size
ep_7 <- fread("epitope_test_7aa_preds.csv", header = T)
ep_7$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_7$class)))
ep_7$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_7$pred)))

ep_9 <- fread("epitope_test_9aa_preds.csv", header = T)
ep_9$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_9$class)))
ep_9$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_9$pred)))

ep_11 <- fread("epitope_test_11aa_preds.csv", header = T)
ep_11$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_11$class)))
ep_11$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_11$pred)))

ep_13 <- fread("epitope_test_13aa_preds.csv", header = T)
ep_13$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_13$class)))
ep_13$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_13$pred)))

ep_15 <- fread("epitope_test_15aa_preds.csv", header = T)
ep_15$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_15$class)))
ep_15$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_15$pred)))

ep_17 <- fread("epitope_test_17aa_preds.csv", header = T)
ep_17$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_17$class)))
ep_17$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_17$pred)))

ep_19 <- fread("epitope_test_19aa_preds.csv", header = T)
ep_19$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_19$class)))
ep_19$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_19$pred)))

ep_21 <- fread("epitope_test_21aa_preds.csv", header = T)
ep_21$class <- as.numeric(gsub("[)]", "", gsub("tensor[(]", "", ep_21$class)))
ep_21$pred <- as.numeric(gsub(",.*", "", gsub("tensor[(]", "", ep_21$pred)))

# load in the ML model 
rf_ep_17 <- fread("17aa_chemistry_rf_test_probabilities.csv", header = T)

# create ROC curves
ep_7_roc <- roc(ep_7$class~ep_7$pred)
ep_9_roc <- roc(ep_9$class~ep_9$pred)
ep_11_roc <- roc(ep_11$class~ep_11$pred)
ep_13_roc <- roc(ep_13$class~ep_13$pred)
ep_15_roc <- roc(ep_15$class~ep_15$pred)
ep_17_roc <- roc(ep_17$class~ep_17$pred)
ep_19_roc <- roc(ep_19$class~ep_19$pred)
ep_21_roc <- roc(ep_21$class~ep_21$pred)

auc_vals <- c(ep_7_roc$auc[1], ep_9_roc$auc[1], ep_11_roc$auc[1], ep_13_roc$auc[1], ep_15_roc$auc[1], ep_17_roc$auc[1], ep_19_roc$auc[1], ep_21_roc$auc[1])
window_size <- seq(7,21,2)
auc_df <- as.data.frame(cbind(window_size, auc_vals))
auc_plot <- ggplot(auc_df, aes(x=window_size, y=auc_vals)) + geom_point() + geom_line() + ylim(.8,1) + scale_x_continuous(breaks=seq(7,21,2))
auc_plot + xlab("Input Window Size (aa)") + ylab("Test Set Performance (AUC)") +theme_bw()+theme(text = element_text(size=16))

rf_ep_17_roc <- roc(rf_ep_17$class~rf_ep_17$probabilty)

# perform roc test comparisons across all DL window sizes
roc_tests <- list(roc.test(ep_7_roc, ep_9_roc, auc=TRUE),
                  roc.test(ep_7_roc, ep_11_roc, auc=TRUE),
                  roc.test(ep_7_roc, ep_13_roc, auc=TRUE),
                  roc.test(ep_7_roc, ep_15_roc, auc=TRUE),
                  roc.test(ep_7_roc, ep_17_roc, auc=TRUE),
                  roc.test(ep_7_roc, ep_19_roc, auc=TRUE),
                  roc.test(ep_7_roc, ep_21_roc, auc=TRUE),
                  roc.test(ep_9_roc, ep_11_roc, auc=TRUE),
                  roc.test(ep_9_roc, ep_13_roc, auc=TRUE),
                  roc.test(ep_9_roc, ep_15_roc, auc=TRUE),
                  roc.test(ep_9_roc, ep_17_roc, auc=TRUE),
                  roc.test(ep_9_roc, ep_19_roc, auc=TRUE),
                  roc.test(ep_9_roc, ep_21_roc, auc=TRUE),
                  roc.test(ep_11_roc, ep_13_roc, auc=TRUE),
                  roc.test(ep_11_roc, ep_15_roc, auc=TRUE),
                  roc.test(ep_11_roc, ep_17_roc, auc=TRUE),
                  roc.test(ep_11_roc, ep_19_roc, auc=TRUE),
                  roc.test(ep_11_roc, ep_21_roc, auc=TRUE),
                  roc.test(ep_13_roc, ep_15_roc, auc=TRUE),
                  roc.test(ep_13_roc, ep_17_roc, auc=TRUE),
                  roc.test(ep_13_roc, ep_19_roc, auc=TRUE),
                  roc.test(ep_13_roc, ep_21_roc, auc=TRUE),
                  roc.test(ep_15_roc, ep_17_roc, auc=TRUE),
                  roc.test(ep_15_roc, ep_19_roc, auc=TRUE),
                  roc.test(ep_15_roc, ep_21_roc, auc=TRUE),
                  roc.test(ep_17_roc, ep_19_roc, auc=TRUE),
                  roc.test(ep_17_roc, ep_21_roc, auc=TRUE),
                  roc.test(ep_19_roc, ep_21_roc, auc=TRUE))

# adjust for multiple comparisons
pvals <- unlist(lapply(roc_tests, function(x){x$p.value}))
pvals_adj <- p.adjust(pvals, method = "BH")
pvals_adj

# compare ML to DL
roc.test(ep_17_roc, rf_ep_17_roc)


# load cross performance on iv data
cross_7 <- fread("ep_mod_cross_performance_7aa.csv", header = T)
cross_17 <- fread("ep_mod_cross_performance_17aa.csv", header = T)

# overall cross
cross_7_roc <- roc(cross_7$labels~cross_7$preds)
cross_17_roc <- roc(cross_17$labels~cross_17$preds)
roc.test(cross_7_roc, cross_17_roc)

cross_7_C <- subset(cross_7, cross_7$proteasome == "C")
cross_17_C <- subset(cross_17, cross_17$proteasome == "C")

cross_7C_roc <- roc(cross_7_C$labels~cross_7_C$preds)
cross_17C_roc <- roc(cross_17_C$labels~cross_17_C$preds)
roc.test(cross_7C_roc, cross_17C_roc)


cross_7_I <- subset(cross_7, cross_7$proteasome == "I")
cross_17_I <- subset(cross_17, cross_17$proteasome == "I")

cross_7I_roc <- roc(cross_7_I$labels~cross_7_I$preds)
cross_17I_roc <- roc(cross_17_I$labels~cross_17_I$preds)
roc.test(cross_7I_roc, cross_17I_roc)

