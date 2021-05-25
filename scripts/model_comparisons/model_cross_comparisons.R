#!usr/bin/env R
# model_cross_comparisons.R
# 
# For issues contact Ben Weeder (weeder@ohsu.edu)
# 
# This script loads trained models across various data types and compares 
# performance on other data types

# set directory where data is stored
setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/validation_results/pepsickle")

# load necessary libraries
library(data.table)
library(ggplot2)
library(pROC)
library(caret)

# cross performance on epitope data
constit_cross <- fread("./pepsickle_constitutive_20S_cross_val_summary.txt")
constit_cross <- constit_cross[order(constit_cross$ID),]
immuno_cross <- fread("./pepsickle_immuno_20S_cross_val_summary.txt")
immuno_cross <- immuno_cross[order(immuno_cross$ID),]

# cross performance on in-vitro data
epi_c_cross <- fread("./pepsickle_epitope_cross_constit_val_summary.txt")
epi_c_cross <- epi_c_cross[order(epi_c_cross$ID),]
epi_i_cross <- fread("./pepsickle_epitope_cross_immuno_val_summary.txt")
epi_i_cross <- epi_i_cross[order(epi_i_cross$ID),]

# in-vitro to in-vitro cross performance
c_i_cross <- fread("./pepsickle_constitutive_20S_cross_val_immuno_summary.txt")
roc(c_i_cross$true_label, c_i_cross$cleavage_prob)
i_c_cross <- fread("./pepsickle_immuno_20S_cross_val_constit_summary.txt")
roc(i_c_cross$true_label, i_c_cross$cleavage_prob)


# performance metrics on epitope data
constit_cross_performance <- caret::confusionMatrix(as.factor(as.numeric(constit_cross$cleavage_pred)), 
                                                as.factor(constit_cross$true_label), mode = "sens_spec", positive="1")
constit_cross_performance
constit_cross_roc <- roc(constit_cross$true_label~constit_cross$cleavage_prob)

immuno_cross_performance <- caret::confusionMatrix(as.factor(as.numeric(immuno_cross$cleavage_pred)), 
                                                    as.factor(immuno_cross$true_label), mode = "sens_spec", positive="1")
immuno_cross_performance
immuno_cross_roc <- roc(immuno_cross$true_label~immuno_cross$cleavage_prob)

# compare i to c cross performance
ci_roc <- roc(c_i_cross$true_label~c_i_cross$cleavage_prob)
ic_roc <- roc(i_c_cross$true_label~i_c_cross$cleavage_prob)

ci_roc
ic_roc
