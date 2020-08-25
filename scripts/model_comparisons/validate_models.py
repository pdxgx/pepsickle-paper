#!/usr/bin/env python3
"""
validate_models.py

For issues contact Ben Weeder (weeder@ohsu.edu)

"""

from sequence_featurization_tools import *
import re
import pickle
import torch
import torch.nn as nn
import pandas as pd
import torch.optim as optim
from sklearn import metrics
import torch.nn.functional as F
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-i", "--in-dir",
                  help="directory of pickled validation data")
parser.add_option("-m", "--model-weights",
                  help="pickled dictionary of trained model weights")
parser.add_option("--human-only", action="store_true",
                  help="flags export files with human only annotation")

(options, args) = parser.parse_args()


class DigestionSeqNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.drop = nn.Dropout(p=0.2)
        self.input = nn.Linear(262, 136)
        self.bn1 = nn.BatchNorm1d(136)
        self.fc1 = nn.Linear(136, 68)
        self.bn2 = nn.BatchNorm1d(68)
        self.fc2 = nn.Linear(68, 34)
        self.bn3 = nn.BatchNorm1d(34)
        self.out = nn.Linear(34, 2)

    def forward(self, x, c_prot, i_prot):
        # make sure input tensor is flattened

        x = x.reshape(x.shape[0], -1)
        x = torch.cat((x, c_prot.reshape(c_prot.shape[0], -1)), 1)
        x = torch.cat((x, i_prot.reshape(i_prot.shape[0], -1)), 1)

        x = self.drop(F.relu(self.bn1(self.input(x))))
        x = self.drop(F.relu(self.bn2(self.fc1(x))))
        x = self.drop(F.relu(self.bn3(self.fc2(x))))
        x = F.log_softmax(self.out(x), dim=1)

        return x


class DigestionMotifNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.drop = nn.Dropout(p=.2)
        self.conv = nn.Conv1d(4, 4, 3, groups=4)
        # self.fc1 = nn.Linear(78, 38)
        self.fc1 = nn.Linear(46, 38)
        self.bn1 = nn.BatchNorm1d(38)
        self.fc2 = nn.Linear(38, 20)
        self.bn2 = nn.BatchNorm1d(20)
        self.out = nn.Linear(20, 2)

    def forward(self, x, c_prot, i_prot):
        # perform convolution prior to flattening
        x = x.transpose(1, 2)
        x = self.conv(x)

        # make sure input tensor is flattened
        x = x.reshape(x.shape[0], -1)
        x = torch.cat((x, c_prot.reshape(c_prot.shape[0], -1)), 1)
        x = torch.cat((x, i_prot.reshape(i_prot.shape[0], -1)), 1)

        x = self.drop(F.relu(self.bn1(self.fc1(x))))
        x = self.drop(F.relu(self.bn2(self.fc2(x))))
        x = F.log_softmax(self.out(x), dim=1)

        return x


class EpitopeSeqNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.drop = nn.Dropout(p=0.2)
        self.input = nn.Linear(260, 136)
        self.bn1 = nn.BatchNorm1d(136)
        self.fc1 = nn.Linear(136, 68)
        self.bn2 = nn.BatchNorm1d(68)
        self.fc2 = nn.Linear(68, 34)
        self.bn3 = nn.BatchNorm1d(34)
        self.out = nn.Linear(34, 2)

    def forward(self, x):
        # make sure input tensor is flattened
        x = x.reshape(x.shape[0], -1)

        x = self.drop(F.relu(self.bn1(self.input(x))))
        x = self.drop(F.relu(self.bn2(self.fc1(x))))
        x = self.drop(F.relu(self.bn3(self.fc2(x))))
        x = F.log_softmax(self.out(x), dim=1)

        return x


class EpitopeMotifNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.drop = nn.Dropout(p=0.2)
        self.conv = nn.Conv1d(4, 4, 3, groups=4)
        self.fc1 = nn.Linear(44, 38)
        self.bn1 = nn.BatchNorm1d(38)
        self.fc2 = nn.Linear(38, 20)
        self.bn2 = nn.BatchNorm1d(20)
        self.out = nn.Linear(20, 2)

    def forward(self, x):
        # perform convolution prior to flattening
        x = x.transpose(1, 2)
        x = self.conv(x)

        # make sure input tensor is flattened
        x = x.reshape(x.shape[0], -1)

        x = self.drop(F.relu(self.bn1(self.fc1(x))))
        x = self.drop(F.relu(self.bn2(self.fc2(x))))
        x = F.log_softmax(self.out(x), dim=1)

        return x


# initialize models
dtype = torch.FloatTensor
model_weights = pickle.load(open(options.model_weights, "rb"))

epitope_motif_model = EpitopeMotifNet().eval()
epitope_sequence_model = EpitopeSeqNet().eval()

digestion_20S_motif_model = DigestionMotifNet().eval()
digestion_20S_sequence_model = DigestionSeqNet().eval()

digestion_26S_motif_model = DigestionMotifNet().eval()
digestion_26S_sequence_model = DigestionSeqNet().eval()

if options.human_only:
    epitope_motif_model.load_state_dict(model_weights['human_epitope_motif_mod'])
    epitope_sequence_model.load_state_dict(model_weights['human_epitope_sequence_mod'])

    digestion_20S_motif_model.load_state_dict(model_weights["human_20S_digestion_motif_mod"])
    digestion_20S_sequence_model.load_state_dict((model_weights["human_20S_digestion_sequence_mod"]))

    digestion_26S_motif_model.load_state_dict((model_weights["human_26S_digestion_motif_mod"]))
    digestion_26S_sequence_model.load_state_dict((model_weights["human_26S_digestion_sequence_mod"]))

else:
    epitope_motif_model.load_state_dict(model_weights['all_mammal_epitope_motif_mod'])
    epitope_sequence_model.load_state_dict(model_weights['all_mammal_epitope_sequence_mod'])

    digestion_20S_motif_model.load_state_dict(model_weights["all_mammal_20S_digestion_motif_mod"])
    digestion_20S_sequence_model.load_state_dict((model_weights["all_mammal_20S_digestion_sequence_mod"]))

    digestion_26S_motif_model.load_state_dict((model_weights["all_mammal_26S_digestion_motif_mod"]))
    digestion_26S_sequence_model.load_state_dict((model_weights["all_mammal_26S_digestion_sequence_mod"]))


## load in and run on epitope data
epitope_val_dict = pickle.load(
    open(options.in_dir + "/epitope_val_filtered.pickle", "rb"))
epitope_positives = list(epitope_val_dict['positives'].keys())
epitope_positive_features = torch.from_numpy(
    generate_feature_array(epitope_positives)).type(dtype)
epitope_negatives = list(epitope_val_dict['negatives'].keys())
epitope_negative_features = torch.from_numpy(
    generate_feature_array(epitope_negatives)).type(dtype)

epitope_pos_preds = torch.exp(
    (epitope_sequence_model(epitope_positive_features[:, :, :20])[:, 1] +
     epitope_motif_model(epitope_positive_features[:, :, 22:])[:, 1]) / 2)

epitope_neg_preds = torch.exp(
    (epitope_sequence_model(epitope_negative_features[:, :, :20])[:, 1] +
     epitope_motif_model(epitope_negative_features[:, :, 22:])[:, 1]) / 2)

epitope_preds = [float(x) for x in epitope_pos_preds] + \
                [float(x) for x in epitope_neg_preds]
epitope_true_labels = [1] * len(epitope_pos_preds) + \
                      [0] * len(epitope_neg_preds)

print("Epitope Validation AUC: ", metrics.roc_auc_score(epitope_true_labels,
                                                        epitope_preds))
print("Epitope Classification Report:\n", metrics.classification_report(
    epitope_true_labels, [x > 0.5 for x in epitope_preds]))
print("\n")


## load in and run on constitutive proteasome data
constitutive_20S_val_dict = pickle.load(
    open(options.in_dir +
         "/20S_digestion_constitutive_validation_filtered.pickle", "rb"))
constitutive_20S_positives = list(constitutive_20S_val_dict['positives'].keys())
constitutive_20S_positive_features = torch.from_numpy(
    generate_feature_array(constitutive_20S_positives)).type(dtype)
constitutive_20S_negatives = list(constitutive_20S_val_dict['negatives'].keys())
constitutive_20S_negative_features = torch.from_numpy(
    generate_feature_array(constitutive_20S_negatives)).type(dtype)

constitutive_20S_pos_preds = torch.exp(
    (digestion_20S_sequence_model(constitutive_20S_positive_features[:, :, :20],
                                  torch.tensor([1] * len(constitutive_20S_positives)).type(dtype),
                                  torch.tensor([0] * len(constitutive_20S_positives)).type(dtype))[:, 1] +
     digestion_20S_motif_model(constitutive_20S_positive_features[:, :, 22:],
                               torch.tensor([1] * len(constitutive_20S_positives)).type(dtype),
                               torch.tensor([0] * len(constitutive_20S_positives)).type(dtype))[:, 1]) / 2)

constitutive_20S_neg_preds = torch.exp(
    (digestion_20S_sequence_model(constitutive_20S_negative_features[:, :, :20],
                                  torch.tensor([1] * len(constitutive_20S_negatives)).type(dtype),
                                  torch.tensor([0] * len(constitutive_20S_negatives)).type(dtype))[:, 1] +
     digestion_20S_motif_model(constitutive_20S_negative_features[:, :, 22:],
                               torch.tensor([1] * len(constitutive_20S_negatives)).type(dtype),
                               torch.tensor([0] * len(constitutive_20S_negatives)).type(dtype))[:, 1]) / 2)

constitutive_20S_preds = [float(x) for x in constitutive_20S_pos_preds] + \
                [float(x) for x in constitutive_20S_neg_preds]
constitutive_20S_true_labels = [1] * len(constitutive_20S_pos_preds) + \
                      [0] * len(constitutive_20S_neg_preds)

print("Constitutive 20S Validation AUC: ", metrics.roc_auc_score(constitutive_20S_true_labels,
                                                        constitutive_20S_preds))
print("constitutive 20S Classification Report:\n", metrics.classification_report(
    constitutive_20S_true_labels, [x > 0.5 for x in constitutive_20S_preds]))
print("\n")


## load in and run on immuno proteasome data
immuno_20S_val_dict = pickle.load(
    open(options.in_dir +
         "/20S_digestion_immuno_validation_filtered.pickle", "rb"))
immuno_20S_positives = list(immuno_20S_val_dict['positives'].keys())
immuno_20S_positive_features = torch.from_numpy(
    generate_feature_array(immuno_20S_positives)).type(dtype)
immuno_20S_negatives = list(immuno_20S_val_dict['negatives'].keys())
immuno_20S_negative_features = torch.from_numpy(
    generate_feature_array(immuno_20S_negatives)).type(dtype)

immuno_20S_pos_preds = torch.exp(
    (digestion_20S_sequence_model(immuno_20S_positive_features[:, :, :20],
                                  torch.tensor([0] * len(immuno_20S_positives)).type(dtype),
                                  torch.tensor([1] * len(immuno_20S_positives)).type(dtype))[:, 1] +
     digestion_20S_motif_model(immuno_20S_positive_features[:, :, 22:],
                               torch.tensor([0] * len(immuno_20S_positives)).type(dtype),
                               torch.tensor([1] * len(immuno_20S_positives)).type(dtype))[:, 1]) / 2)

immuno_20S_neg_preds = torch.exp(
    (digestion_20S_sequence_model(immuno_20S_negative_features[:, :, :20],
                                  torch.tensor([0] * len(immuno_20S_negatives)).type(dtype),
                                  torch.tensor([1] * len(immuno_20S_negatives)).type(dtype))[:, 1] +
     digestion_20S_motif_model(immuno_20S_negative_features[:, :, 22:],
                               torch.tensor([0] * len(immuno_20S_negatives)).type(dtype),
                               torch.tensor([1] * len(immuno_20S_negatives)).type(dtype))[:, 1]) / 2)

immuno_20S_preds = [float(x) for x in immuno_20S_pos_preds] + \
                [float(x) for x in immuno_20S_neg_preds]
immuno_20S_true_labels = [1] * len(immuno_20S_pos_preds) + \
                      [0] * len(immuno_20S_neg_preds)

print("immunoproteasome 20S Validation AUC: ", metrics.roc_auc_score(immuno_20S_true_labels,
                                                        immuno_20S_preds))
print("immunoproteasome 20S Classification Report:\n", metrics.classification_report(
    immuno_20S_true_labels, [x > 0.5 for x in immuno_20S_preds]))
print("\n")


# assess 26S models
digestion_26S_dat = pd.read_csv(options.in_dir +
                                "/all_mammal_26S_val_data.csv")

constit_26S_dat = digestion_26S_dat.loc[(digestion_26S_dat['c_prot'] == 1) &
                                        (digestion_26S_dat['i_prot'] == 0)]

constit_26S_positives = constit_26S_dat[constit_26S_dat['cleaved'] == 1][
    'window']
constit_26S_negatives = constit_26S_dat[constit_26S_dat['cleaved'] == 0][
    'window'].sample(constit_26S_positives.shape[0])

constitutive_26S_positive_features = torch.from_numpy(
    generate_feature_array(constit_26S_positives)).type(dtype)

constitutive_26S_negative_features = torch.from_numpy(
    generate_feature_array(constit_26S_negatives)).type(dtype)


constitutive_26S_pos_preds = torch.exp(
    (digestion_26S_sequence_model(constitutive_26S_positive_features[:, :, :20],
                                  torch.tensor([1] * len(constit_26S_positives)).type(dtype),
                                  torch.tensor([0] * len(constit_26S_positives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(constitutive_26S_positive_features[:, :, 22:],
                               torch.tensor([1] * len(constit_26S_positives)).type(dtype),
                               torch.tensor([0] * len(constit_26S_positives)).type(dtype))[:, 1]) / 2)

constitutive_26S_neg_preds = torch.exp(
    (digestion_26S_sequence_model(constitutive_26S_negative_features[:, :, :20],
                                  torch.tensor([1] * len(constit_26S_negatives)).type(dtype),
                                  torch.tensor([0] * len(constit_26S_negatives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(constitutive_26S_negative_features[:, :, 22:],
                               torch.tensor([1] * len(constit_26S_negatives)).type(dtype),
                               torch.tensor([0] * len(constit_26S_negatives)).type(dtype))[:, 1]) / 2)

constitutive_26S_preds = [float(x) for x in constitutive_26S_pos_preds] + \
                [float(x) for x in constitutive_26S_neg_preds]
constitutive_26S_true_labels = [1] * len(constitutive_26S_pos_preds) + \
                      [0] * len(constitutive_26S_neg_preds)

print("Constitutive 26S Validation AUC: ", metrics.roc_auc_score(
    constitutive_26S_true_labels, constitutive_26S_preds))
print("constitutive 26S Classification Report:\n", metrics.classification_report(
    constitutive_26S_true_labels, [x > 0.5 for x in constitutive_26S_preds]))
print("\n")


# immuno
immuno_26S_dat = digestion_26S_dat.loc[(digestion_26S_dat['c_prot'] == 0) &
                                       (digestion_26S_dat['i_prot'] == 1)]

immuno_26S_positives = immuno_26S_dat[immuno_26S_dat['cleaved'] == 1][
    'window']
immuno_26S_negatives = immuno_26S_dat[immuno_26S_dat['cleaved'] == 0][
    'window'].sample(immuno_26S_positives.shape[0])


immuno_26S_positive_features = torch.from_numpy(
    generate_feature_array(immuno_26S_positives)).type(dtype)

immuno_26S_negative_features = torch.from_numpy(
    generate_feature_array(immuno_26S_negatives)).type(dtype)


immuno_26S_pos_preds = torch.exp(
    (digestion_26S_sequence_model(immuno_26S_positive_features[:, :, :20],
                                  torch.tensor([0] * len(immuno_26S_positives)).type(dtype),
                                  torch.tensor([1] * len(immuno_26S_positives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(immuno_26S_positive_features[:, :, 22:],
                               torch.tensor([0] * len(immuno_26S_positives)).type(dtype),
                               torch.tensor([1] * len(immuno_26S_positives)).type(dtype))[:, 1]) / 2)

immuno_26S_neg_preds = torch.exp(
    (digestion_26S_sequence_model(immuno_26S_negative_features[:, :, :20],
                                  torch.tensor([0] * len(immuno_26S_negatives)).type(dtype),
                                  torch.tensor([1] * len(immuno_26S_negatives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(immuno_26S_negative_features[:, :, 22:],
                               torch.tensor([0] * len(immuno_26S_negatives)).type(dtype),
                               torch.tensor([1] * len(immuno_26S_negatives)).type(dtype))[:, 1]) / 2)

immuno_26S_preds = [float(x) for x in immuno_26S_pos_preds] + \
                [float(x) for x in immuno_26S_neg_preds]
immuno_26S_true_labels = [1] * len(immuno_26S_pos_preds) + \
                      [0] * len(immuno_26S_neg_preds)

print("Immunoproteasome 26S Validation AUC: ", metrics.roc_auc_score(
    immuno_26S_true_labels, immuno_26S_preds))
print("Immunoproteasome 26S Classification Report:\n", metrics.classification_report(
    immuno_26S_true_labels, [x > 0.5 for x in immuno_26S_preds]))
print("\n")


# 26S model applied to constitutive 20S data
constitutive_20S_26S_cross_pos_preds = torch.exp(
    (digestion_26S_sequence_model(constitutive_20S_positive_features[:, :, :20],
                                  torch.tensor([1] * len(constitutive_20S_positives)).type(dtype),
                                  torch.tensor([0] * len(constitutive_20S_positives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(constitutive_20S_positive_features[:, :, 22:],
                               torch.tensor([1] * len(constitutive_20S_positives)).type(dtype),
                               torch.tensor([0] * len(constitutive_20S_positives)).type(dtype))[:, 1]) / 2)

constitutive_20S_26S_cross_neg_preds = torch.exp(
    (digestion_26S_sequence_model(constitutive_20S_negative_features[:, :, :20],
                                  torch.tensor([1] * len(constitutive_20S_negatives)).type(dtype),
                                  torch.tensor([0] * len(constitutive_20S_negatives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(constitutive_20S_negative_features[:, :, 22:],
                               torch.tensor([1] * len(constitutive_20S_negatives)).type(dtype),
                               torch.tensor([0] * len(constitutive_20S_negatives)).type(dtype))[:, 1]) / 2)

constitutive_20S_26S_cross_preds = [float(x) for x in constitutive_20S_26S_cross_pos_preds] + \
                [float(x) for x in constitutive_20S_26S_cross_neg_preds]

print("Constitutive 20S by 26S model Validation AUC: ",
      metrics.roc_auc_score(constitutive_20S_true_labels,
                            constitutive_20S_26S_cross_preds))
print("constitutive 20S by 26S model Classification Report:\n",
      metrics.classification_report(constitutive_20S_true_labels,
                                    [x > 0.5 for x in constitutive_20S_26S_cross_preds]))
print("\n")


# 26S model on 20S immunoproteasome data
immuno_20S_26S_cross_pos_preds = torch.exp(
    (digestion_26S_sequence_model(immuno_20S_positive_features[:, :, :20],
                                  torch.tensor([0] * len(immuno_20S_positives)).type(dtype),
                                  torch.tensor([1] * len(immuno_20S_positives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(immuno_20S_positive_features[:, :, 22:],
                               torch.tensor([0] * len(immuno_20S_positives)).type(dtype),
                               torch.tensor([1] * len(immuno_20S_positives)).type(dtype))[:, 1]) / 2)

immuno_20S_26S_cross_neg_preds = torch.exp(
    (digestion_26S_sequence_model(immuno_20S_negative_features[:, :, :20],
                                  torch.tensor([0] * len(immuno_20S_negatives)).type(dtype),
                                  torch.tensor([1] * len(immuno_20S_negatives)).type(dtype))[:, 1] +
     digestion_26S_motif_model(immuno_20S_negative_features[:, :, 22:],
                               torch.tensor([0] * len(immuno_20S_negatives)).type(dtype),
                               torch.tensor([1] * len(immuno_20S_negatives)).type(dtype))[:, 1]) / 2)

immuno_20S_26S_cross_preds = [float(x) for x in immuno_20S_26S_cross_pos_preds] + \
                [float(x) for x in immuno_20S_26S_cross_neg_preds]

print("immunoproteasome 20S by 26S model Validation AUC: ",
      metrics.roc_auc_score(immuno_20S_true_labels, immuno_20S_26S_cross_preds))
print("immunoproteasome 20S by 26S model Classification Report:\n",
      metrics.classification_report(immuno_20S_true_labels,
                                    [x > 0.5 for x in immuno_20S_26S_cross_preds]))
print("\n")
