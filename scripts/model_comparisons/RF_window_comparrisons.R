#!usr/bin/env R

# ML_window_comparisons.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads windows of various sizes for predicted cleavage sites and compares test performances

# Load in required packages
# setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/in-vitro_mod_window_preds")
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/output/in-vitro_mod_window_preds")

library(pROC)
library(data.table)

# load in 7aa window predictions (ML)
iv_7aa <- fread("7aa_chemistry_rf_test_probabilities.csv", header = T)
iv_7aa_C <- subset(iv_7aa, iv_7aa$proteasome == "C")
iv_7aa_I <- subset(iv_7aa, iv_7aa$proteasome == "I")

# load in 21aa window predictions (ML)
iv_21aa <- fread("21aa_chemistry_rf_test_probabilities.csv", header = T)
iv_21aa_C <- subset(iv_21aa, iv_21aa$proteasome == "C")
iv_21aa_I <- subset(iv_21aa, iv_21aa$proteasome == "I")

# load in 7aa window dredictions (DL)
iv_net_7aa <- fread("in_vitro_net_test_7aa_preds.csv")
colnames(iv_net_7aa)[4] <- "probability"

# look at in-vitro constitutive ML performance based on size
iv_7aa_Croc <- roc(iv_7aa_C$class~iv_7aa_C$probabilty)
iv_21aa_Croc <- roc(iv_21aa_C$class~iv_21aa_C$probabilty)
roc.test(iv_7aa_Croc, iv_21aa_Croc, auc=T)

# look at in-vitro immuno ML performance based on size
iv_7aa_Iroc <- roc(iv_7aa_I$class~iv_7aa_I$probabilty)
iv_21aa_Iroc <- roc(iv_21aa_I$class~iv_21aa_I$probabilty)
roc.test(iv_7aa_Iroc, iv_21aa_Iroc, auc=T)

# compare overall performance based on size
overall_7aa_roc <- roc(iv_7aa$class~iv_7aa$probabilty)
overall_21aa_roc <- roc(iv_21aa$class~iv_21aa$probabilty)
roc.test(overall_7aa_roc, overall_21aa_roc)

# compare 7aa overall for DL and ML
overall_net_7aa_roc <- roc(iv_net_7aa$class~iv_net_7aa$probability)
roc.test(overall_7aa_roc, overall_net_7aa_roc)

