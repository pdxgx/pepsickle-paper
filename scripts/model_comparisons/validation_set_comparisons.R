#!usr/bin/env R

# validation_set_comparisons.R
#
# For issues contact Ben Weeder (weeder@ohsu.edu)
#
# This script compares validation performance across tools

# setwd("~/PycharmProjects/pepsickle-paper/data/validation_data/validation_results")
setwd("./pepsickle-paper/data/validation_data/validation_results")

library(data.table)
library(ggplot2)
library(pROC)
library(caret)


# tmp
tmp_dat <- fread("~/PycharmProjects/pepsickle-paper/data/model_weights/final_test_preds.csv")
ggplot(tmp_dat, aes(x=cleave_p, fill=as.factor(true_label), group=as.factor(true_label)))+geom_density(alpha = .5)+xlim(0,1)+geom_vline(xintercept = .361)

## Load in netchop data
# epitope
netchop_epitope_val <- fread("./netchop/parsed_netchop_epitope_val_21aa_cleavage_preds.txt")
netchop_epitope_val <- netchop_epitope_val[order(netchop_epitope_val$ident),]
# 20S immunoproteasome
netchop_immuno_val <- fread("./netchop/parsed_netchop_20S_immuno_digestion_val_21aa_cleavage_preds.txt")
netchop_immuno_val <- netchop_immuno_val[order(netchop_immuno_val$ident),]
# 20S constitutive proteasome
netchop_constit_val <- fread("./netchop/parsed_netchop_20S_constitutive_digestion_val_21aa_cleavage_preds.txt")
netchop_constit_val <- netchop_constit_val[order(netchop_constit_val$ident),]

## Load in PCPS data
# epitope
pcps_epitope_val <- fread("./PCPS/PCPS_epitope_pred_summary.csv")
pcps_epitope_val$ID <- tolower(pcps_epitope_val$ID) # edit ID's to match others
pcps_epitope_val <- pcps_epitope_val[order(pcps_epitope_val$ID),]

# constitutive data
pcps_constit_val <- fread("./PCPS/PCPS_20S_constit_pred_summary.csv")
pcps_constit_val$ID <- tolower(pcps_constit_val$ID)
pcps_constit_val <- pcps_constit_val[order(pcps_constit_val$ID),]

# immuno data
pcps_immuno_val <- fread("./PCPS/PCPS_20S_immuno_summary.csv")
pcps_immuno_val$ID <- tolower(pcps_immuno_val$ID)
pcps_immuno_val <- pcps_immuno_val[order(pcps_immuno_val$ID),]

## Load in pepsickle data
# epitope
pepsickle_epitope_val <- fread("pepsickle/pepsickle_epitope_val_summary.txt")
pepsickle_epitope_val <- pepsickle_epitope_val[order(pepsickle_epitope_val$ID),]
# constit
pepsickle_constit_val <- fread("pepsickle/pepsickle_constitutive_20S_val_summary.txt")
pepsickle_constit_val <- pepsickle_constit_val[order(pepsickle_constit_val$ID),]
#immuno
pepsickle_immuno_val <- fread("pepsickle/pepsickle_immuno_20S_val_summary.txt")
pepsickle_immuno_val <- pepsickle_immuno_val[order(pepsickle_immuno_val$ID),]

## add in PCleavage here
pcleavage_epitope_val <- fread("Pcleavage/PCleavage_21aa_epitope_results.txt")
colnames(pcleavage_epitope_val) <- c("ID", "sequence", "true_label", "score", "cleavage_pred")
pcleavage_epitope_val <- pcleavage_epitope_val[order(pcleavage_epitope_val$ID),]

## NOTE: Predicted with epitope model due to size restrictions
pcleavage_constit_val <- fread("Pcleavage/pcleavage_20S_constitutive_val_preds.txt")
colnames(pcleavage_constit_val) <- c("ID", "sequence", "true_label", "score", "cleavage_pred")
pcleavage_constit_val <- pcleavage_constit_val[order(pcleavage_constit_val$ID),]

pcleavage_immuno_val <- fread("Pcleavage/pcleavage_20S_immuno_val_preds.txt")
colnames(pcleavage_immuno_val) <- c("ID", "sequence", "true_label", "score", "cleavage_pred")
pcleavage_immuno_val <- pcleavage_immuno_val[order(pcleavage_immuno_val$ID),]


## Compare epitope roc curves
netchop_epitope_roc <- roc(netchop_epitope_val$true_class~netchop_epitope_val$prob)
pepsickle_epitope_roc <- roc(pepsickle_epitope_val$true_label~pepsickle_epitope_val$cleavage_prob)
pcleavage_epitope_roc <- roc(pcleavage_epitope_val$true_label~pcleavage_epitope_val$score)

# pick best of both models given
pcps_epitope_roc_c <- roc(pcps_epitope_val$true_label~pcps_epitope_val$constit_prob) # better
pcps_epitope_roc_i <- roc(pcps_epitope_val$true_label~pcps_epitope_val$immuno_prob)
roc.test(pcps_epitope_roc_c, pcps_epitope_roc_i)

# compare pepsickle to the other models available
p.adjust(c(roc.test(pepsickle_epitope_roc, netchop_epitope_roc, paired = TRUE)$p.value, 
         roc.test(pepsickle_epitope_roc, pcps_epitope_roc_c, paired = TRUE)$p.value,
         roc.test(pepsickle_epitope_roc, pcleavage_epitope_roc, paired = TRUE)$p.value),
         method = "BH")


pepsickle_performance <- caret::confusionMatrix(as.factor(as.numeric(pepsickle_epitope_val$cleavage_pred)), 
                       as.factor(pepsickle_epitope_val$true_label), mode = "prec_recall", positive="1")
# for sensitivity comparison
caret::confusionMatrix(as.factor(as.numeric(pepsickle_epitope_val$cleavage_pred)), 
                       as.factor(pepsickle_epitope_val$true_label), mode = "sens_spec", positive="1")


netchop_performance <- caret::confusionMatrix(as.factor(netchop_epitope_val$cleave_pred), 
                                                as.factor(netchop_epitope_val$true_class), mode = "prec_recall", positive="1")
pcps_performance <- caret::confusionMatrix(as.factor(as.numeric(pcps_epitope_val$constit_mod_pred)), 
                                           as.factor(pcps_epitope_val$true_label), mode = "prec_recall", positive="1")
pcleavage_performance <- caret::confusionMatrix(as.factor(as.numeric(pcleavage_epitope_val$cleavage_pred)), 
                                           as.factor(pcleavage_epitope_val$true_label), mode = "prec_recall", positive="1")


epitope_out <- rbind(pepsickle_performance$byClass[c("Precision", "Recall", "F1")],
                     netchop_performance$byClass[c("Precision", "Recall", "F1")],
                     pcps_performance$byClass[c("Precision", "Recall", "F1")],
                     pcleavage_performance$byClass[c("Precision", "Recall", "F1")])
rownames(epitope_out) <- c("Pepsickle", "NetChop", "PCPS", "PCleavage")
write.csv(epitope_out, file = "epitope_model_performance_scores.csv")

# plot epitope AUC's
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
epitope_plot <- ggroc(list(pepsickle_epitope_roc, netchop_epitope_roc, pcleavage_epitope_roc, pcps_epitope_roc_c), aes = c("color"), size = 1.2) + 
  theme_classic() + geom_abline(intercept = 1, slope = 1, linetype=2)+
  scale_color_manual(values=c(cbPalette[7], cbPalette[3], cbPalette[2], cbPalette[4]), 
                     labels=c("Pepsickle (0.878)", "Netchop (0.793)", "PCleavage (0.772)", "PCPS (0.761)"))+
  # scale_linetype_manual(values = c("longdash", "4C88C488", "dotdash", "solid"))+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(legend.position = c(.85, .15)) + labs(colour="Model") + 
  # ggtitle("Figure X. Epitope validation performance.") + 
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size=12),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size=12))
epitope_plot
ggsave("~/PycharmProjects/pepsickle-paper/figures/main/epitope_auc.png", epitope_plot, device = "png", width = 8, height = 6, units = "in")

## Compare constit roc curves
netchop_constit_roc <- roc(netchop_constit_val$true_class~netchop_constit_val$prob)
pepsickle_constit_roc <- roc(pepsickle_constit_val$true_label~pepsickle_constit_val$cleavage_prob)
pcleavage_constit_roc <- roc(pcleavage_constit_val$true_label~pcleavage_constit_val$score)

# pick best of both models given
pcps_constit_roc_c <- roc(pcps_constit_val$true_label~pcps_constit_val$constit_prob) # better
pcps_constit_roc_i <- roc(pcps_constit_val$true_label~pcps_constit_val$immuno_prob)
roc.test(pcps_constit_roc_c, pcps_constit_roc_i)
# compare pepsickle to the other models available
roc.test(pepsickle_constit_roc, netchop_constit_roc, paired = TRUE)
roc.test(pepsickle_constit_roc, pcps_constit_roc_c, paired = TRUE)
roc.test(pepsickle_constit_roc, pcleavage_constit_roc, paired = TRUE)
constit_p <- p.adjust(c(roc.test(pepsickle_constit_roc, netchop_constit_roc, paired = TRUE)$p.value,
         roc.test(pepsickle_constit_roc, pcps_constit_roc_c, paired = TRUE)$p.value),
         method = "BH")
round(constit_p, 3) # sig

pepsickle_constit_val$old_prob <- pepsickle_constit_val$cleavage_prob - (.5 - .361)
ggplot(pepsickle_constit_val, aes(x=cleavage_prob, fill=as.factor(true_label), group=as.factor(true_label)))+geom_density(alpha = .5)+xlim(0,1)+geom_vline(xintercept = .361)

pepsickle_c_performance <- caret::confusionMatrix(as.factor(as.numeric(pepsickle_constit_val$cleavage_pred)), 
                                                as.factor(pepsickle_constit_val$true_label), mode = "prec_recall", positive="1")
netchop_c_performance <- caret::confusionMatrix(as.factor(netchop_constit_val$cleave_pred), 
                                              as.factor(netchop_constit_val$true_class), mode = "prec_recall", positive="1")
pcps_c_performance <- caret::confusionMatrix(as.factor(as.numeric(pcps_constit_val$constit_mod_pred)), 
                                           as.factor(pcps_constit_val$true_label), mode = "prec_recall", positive="1")

constit_out <- rbind(pepsickle_c_performance$byClass[c("Precision", "Recall", "F1")],
                     netchop_c_performance$byClass[c("Precision", "Recall", "F1")],
                     pcps_c_performance$byClass[c("Precision", "Recall", "F1")])
rownames(constit_out) <- c("Pepsickle", "NetChop", "PCPS")
write.csv(constit_out, file = "constit_model_performance_scores.csv")

# plot epitope AUC's
constit_plot <- ggroc(list(pepsickle_constit_roc, netchop_constit_roc, pcps_constit_roc_c), aes = c("color"), size = 1.2) + 
  theme_classic() + geom_abline(intercept = 1, slope = 1, linetype=2)+
  scale_color_manual(values=c(cbPalette[7], cbPalette[3], cbPalette[4]), 
                     labels=c("Pepsickle (0.821)", "Netchop (0.801)", "PCPS (0.584)"))+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(legend.position = c(.85, .15)) + labs(colour="Model") + 
  theme(legend.title = element_text(size = 14), legend.text = element_text(size=12), axis.title = element_text(size = 16), axis.text = element_text(size=12))+
  # ggtitle("Figure X. Consitutive Proteasome validation performance") +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size=12),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size=12))
constit_plot
ggsave("~/PycharmProjects/pepsickle-paper/figures/main/constit_auc.png", constit_plot, device = "png", width = 8, height = 6, units = "in")


## Compare immuno roc curves
netchop_immuno_roc <- roc(netchop_immuno_val$true_class~netchop_immuno_val$prob)
pepsickle_immuno_roc <- roc(pepsickle_immuno_val$true_label~pepsickle_immuno_val$cleavage_prob)
pcleavage_immuno_roc <- roc(pcleavage_immuno_val$true_label~pcleavage_immuno_val$score)
# pick best of both models given
pcps_immuno_roc_c <- roc(pcps_immuno_val$true_label~pcps_immuno_val$constit_prob) # better
pcps_immuno_roc_i <- roc(pcps_immuno_val$true_label~pcps_immuno_val$immuno_prob)
roc.test(pcps_immuno_roc_c, pcps_immuno_roc_i)
# compare pepsickle to the other models available
roc.test(pepsickle_immuno_roc, netchop_immuno_roc, paired = TRUE)
roc.test(pepsickle_immuno_roc, pcps_immuno_roc_i, paired = TRUE)
roc.test(pepsickle_immuno_roc, pcleavage_immuno_roc, paired = TRUE)

p.adjust(c(roc.test(pepsickle_immuno_roc, netchop_immuno_roc, paired = TRUE)$p.value,
           roc.test(pepsickle_immuno_roc, pcps_immuno_roc_i, paired = TRUE)$p.value),
         method = "BH")

ggplot(pepsickle_immuno_val, aes(x=cleavage_prob, fill=as.factor(true_label), group=as.factor(true_label)))+geom_density(alpha = .5)+geom_vline(xintercept = .361)

pepsickle_immuno_val$p2 <- pepsickle_immuno_val$cleavage_prob > .361

pepsickle_i_performance <- caret::confusionMatrix(as.factor(as.numeric(pepsickle_immuno_val$cleavage_pred)), 
                                                  as.factor(pepsickle_immuno_val$true_label), mode = "prec_recall", positive="1")
netchop_i_performance <- caret::confusionMatrix(as.factor(netchop_immuno_val$cleave_pred), 
                                                as.factor(netchop_immuno_val$true_class), mode = "prec_recall", positive="1")
pcps_i_performance <- caret::confusionMatrix(as.factor(as.numeric(pcps_immuno_val$constit_mod_pred)), 
                                             as.factor(pcps_immuno_val$true_label), mode = "prec_recall", positive="1")

immuno_out <- rbind(pepsickle_i_performance$byClass[c("Precision", "Recall", "F1")],
                     netchop_i_performance$byClass[c("Precision", "Recall", "F1")],
                     pcps_i_performance$byClass[c("Precision", "Recall", "F1")])
rownames(immuno_out) <- c("Pepsickle", "NetChop", "PCPS")
write.csv(immuno_out, file = "immuno_model_performance_scores.csv")


# plot immuno AUC's
immuno_plot <- ggroc(list(pepsickle_immuno_roc, netchop_immuno_roc, pcps_immuno_roc_i), aes = c("color"), size = 1.2) + 
  theme_classic() + geom_abline(intercept = 1, slope = 1, linetype=2)+
  scale_color_manual(values=c(cbPalette[7], cbPalette[3], cbPalette[4]), 
                     labels=c("Pepsickle (0.789)", "Netchop (0.659)", "PCPS (0.623)"))+
  xlab("Specificity")+
  ylab("Sensitivity")+
  theme(legend.position = c(.85, .15)) + labs(colour="Model") + 
  theme(legend.title = element_text(size = 14), legend.text = element_text(size=12), axis.title = element_text(size = 16), axis.text = element_text(size=12))+
  # ggtitle("Figure X. Immuno Proteasome validation performance")+ 
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size=12),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size=12))
immuno_plot
ggsave("~/PycharmProjects/pepsickle-paper/figures/main/immuno_auc.png", immuno_plot, device = "png", width = 8, height = 6, units = "in")

immuno_plot


epitope_plot <- ggroc(list(pepsickle_epitope_roc, netchop_epitope_roc, pcps_epitope_roc_c), aes = c("color"), size = 1.2) + 
  theme_classic() + geom_abline(intercept = 1, slope = 1, linetype=2)+
  scale_color_manual(values=c("#e62e1e", "#2975f0", "#07a60c"), labels=c("Pepsickle Epitope Model", "Netchop Epitope Model", "PCPS Model"))+
  theme(legend.position = "bottom") + labs(colour="Model") + 
  ggtitle("Figure X. Epitope validation performance") + 
  theme(legend.title = element_text(size = 14), legend.text = element_text(size=12), axis.title = element_text(size = 14))
epitope_plot


all_plot <- ggroc(list(constit_roc, immuno_roc, netchop_constit_roc, netchop_immuno_roc), aes("linetype", "color")) +
  geom_abline(intercept = 1, slope = 1, linetype=2)+
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed"), labels=c("Pepsickle Model (constitutive)", 
                                                                      "Pepsickle Model (immuno)",
                                                                      "NetChop Model (constitutive)", 
                                                                      "NetChop Model (immuno)")) + 
  scale_color_manual(values = c("red", "red", "blue", "blue")) +
  guides(color=FALSE)
all_plot


power.roc.test(pepsickle_immuno_roc, netchop_immuno_roc, sig.level = 0.05)
power.roc.test(pepsickle_immuno_roc, netchop_immuno_roc, power = 0.8)
