#!/usr/bin/env python3
"""
epitope_based_ensemble_net.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script trains neural networks based on both the sequence identity
and physical property motifs of cleavage and non-cleavage examples from epitope
databases. Exports trained model wieghts.
"""
from sequence_featurization_tools import *
from captum.attr import Saliency
import pickle
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn import metrics
import torch.nn.functional as F
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-i", "--in-file",
                  help="pickled dictionary of cleavage windows")
parser.add_option("-o", "--out",
                  help="output directory where results will be exported")
parser.add_option("-w", "--window-size", default=17,
                  help="sets input wintow size for DL models")
parser.add_option("--human-only", action="store_true",
                  help="flags export files with human only annotation")
parser.add_option("--GPU", action="store_true",
                  help="trains models using available GPU")
(options, args) = parser.parse_args()

# set CPU or GPU training
if options.GPU:
    dtype = torch.cuda.FloatTensor
else:
    dtype = torch.FloatTensor

test_holdout_p = .2  # proportion of data held out for testing set
n_epoch = 36

# set seed for consistency
torch.manual_seed(123)


# define model structures
class SeqNet(nn.Module):
    def __init__(self):
        super().__init__()
        # self.in_nodes = 260 # for normal 13aa window
        self.in_nodes = int(options.window_size) * 20
        self.drop = nn.Dropout(p=0.2)
        self.input = nn.Linear(self.in_nodes, 136)
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


class MotifNet(nn.Module):
    def __init__(self):
        super().__init__()
        # self.in_nodes = 44
        self.in_nodes = (int(options.window_size) - 2) * 4
        self.drop = nn.Dropout(p=0.2)
        self.conv = nn.Conv1d(4, 4, 3, groups=4)
        self.fc1 = nn.Linear(self.in_nodes, 38)
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


# initialize networks
sequence_model = SeqNet()
motif_model = MotifNet()

# convert models to cuda if on GPU
if dtype is torch.cuda.FloatTensor:
    sequence_model = sequence_model.cuda()
    motif_model = motif_model.cuda()

#  open pickled dictionary and load in data
handle = open(options.in_file, "rb")
data = pickle.load(handle)
# subset to epitope data for this model
data = data['epitope']

# create list of cleavage windows
pos_windows = []
for key in data['positives'].keys():
    pos_windows.append(key)

# generate features
pos_features = generate_feature_array(pos_windows)
pos_feature_matrix = torch.from_numpy(pos_features)

# create list of non cleavage windows
neg_windows = []
neg_digestion_windows = []
for key in data['negatives'].keys():
    neg_windows.append(key)

# generate features
neg_features = generate_feature_array(neg_windows)
neg_feature_matrix = torch.from_numpy(neg_features)

# define number of training cases based on holdout (unbalanced)
pos_train_k = round((1-test_holdout_p) * pos_feature_matrix.size(0))
neg_train_k = round((1-test_holdout_p) * neg_feature_matrix.size(0))

# permute and split data
pos_perm = torch.randperm(pos_feature_matrix.size(0))
pos_train = pos_feature_matrix[pos_perm[:pos_train_k]]
pos_test = pos_feature_matrix[pos_perm[pos_train_k:]]

neg_perm = torch.randperm(neg_feature_matrix.size(0))
neg_train = neg_feature_matrix[neg_perm[:neg_train_k]]
# use same number of negative test examples as positives for balanced set
neg_test = neg_feature_matrix[neg_perm[neg_train_k:(pos_test.size(0) +
                                                    neg_train_k)]]

# pair training data with labels
pos_train_labeled = []
for i in range(len(pos_train)):
    pos_train_labeled.append([pos_train[i], torch.tensor(1)])
neg_train_labeled = []
for i in range(len(neg_train)):
    neg_train_labeled.append([neg_train[i], torch.tensor(0)])

# combine cleavage and non-cleavage train examples
train_data = pos_train_labeled + neg_train_labeled
train_loader = torch.utils.data.DataLoader(train_data, batch_size=64,
                                           shuffle=True)

# pair test data with labels
pos_test_labeled = []
for i in range(len(pos_test)):
    pos_test_labeled.append([pos_test[i], torch.tensor(1)])
neg_test_labeled = []
for i in range(len(neg_test)):
    neg_test_labeled.append([neg_test[i], torch.tensor(0)])

# combine cleavage and non-cleavage test examples
test_data = pos_test_labeled + neg_test_labeled
# set batch size to give unique set with no re-use
test_loader = torch.utils.data.DataLoader(
    test_data, batch_size=int(len(test_data)/n_epoch + 1), shuffle=True)


m, n, r = pos_features[:, :, 22:].shape
ep_training_features = pd.DataFrame(pos_features[:, :, 22:].reshape(-1, n*r))
ep_training_features['class'] = "pos"
tmp_neg = pd.DataFrame(neg_features[:, :, 22:].reshape(-1, n*r))
tmp_neg['class'] = "neg"
ep_training_features = ep_training_features.append(tmp_neg)
ep_training_features.to_csv("/Users/weeder/PycharmProjects/pepsickle-paper/"
                            "data/validation_data/output/plots/"
                            "epitope_physical_training_windows" +
                            str(options.window_size) + "aa.csv")
print("physical training windows exported")

m, n, r = pos_features[:, :, :20].shape
ep_training_features = pd.DataFrame(pos_features[:, :, :20].reshape(-1, n*r))
ep_training_features['class'] = "pos"
tmp_neg = pd.DataFrame(neg_features[:, :, :20].reshape(-1, n*r))
tmp_neg['class'] = "neg"
ep_training_features = ep_training_features.append(tmp_neg)
ep_training_features.to_csv("/Users/weeder/PycharmProjects/pepsickle-paper/"
                            "data/validation_data/output/plots/"
                            "epitope_sequence_training_windows" +
                            str(options.window_size) + "aa.csv")
print("sequence training windows exported")


# establish training parameters
# inverse weighting for class imbalance in training set
seq_criterion = nn.NLLLoss(
    weight=torch.tensor([1, len(neg_train)/len(pos_train)]).type(dtype))
seq_optimizer = optim.Adam(sequence_model.parameters(), lr=.001)

motif_criterion = nn.NLLLoss(
    weight=torch.tensor([1, len(neg_train)/len(pos_train)]).type(dtype))
motif_optimizer = optim.Adam(motif_model.parameters(), lr=.001)

# initialize tracking of optimal models and train
prev_seq_auc = 0
prev_motif_auc = 0

for epoch in range(n_epoch):
    # reset running loss for each epoch
    seq_running_loss = 0
    motif_running_loss = 0
    # load data, convert if needed
    for dat, labels in train_loader:
        # convert to proper data type
        dat = dat.type(dtype)
        if dtype == torch.cuda.FloatTensor:
            labels = labels.cuda()

        # reset gradients
        seq_optimizer.zero_grad()
        motif_optimizer.zero_grad()

        # generate model predictions
        seq_est = sequence_model(dat[:, :, :20])  # one hot encoded sequences
        # with torch.no_grad():
        #     motif_dat = conv_pre(dat[:, :, 22:].transpose(1, 2))
        motif_est = motif_model(dat[:, :, 22:])  # physical properties (not side chains)

        # calculate loss
        seq_loss = seq_criterion(seq_est, labels)
        motif_loss = motif_criterion(motif_est, labels)

        # back prop loss and step
        seq_loss.backward()
        seq_optimizer.step()
        motif_loss.backward()
        motif_optimizer.step()

        seq_running_loss += seq_loss.item()
        motif_running_loss += motif_loss.item()
    else:
        # output progress
        print(f'Epoch: {epoch + 1}')
        print(f'Sequence Model Running Loss: {seq_running_loss}')
        print(f'Motif Model Running Loss: {motif_running_loss}')

        # test with no grad to speed up
        with torch.no_grad():
            sequence_model.eval()
            motif_model.eval()
            dat, labels = next(iter(test_loader))
            # convert to proper data type
            dat = dat.type(dtype)

            # get est probability of cleavage event
            exp_seq_est = torch.exp(sequence_model(dat[:, :, :20]))[:, 1].cpu()
            # motif_dat = conv_pre(dat[:, :, 22:].transpose(1, 2))
            exp_motif_est = torch.exp(motif_model(dat[:, :, 22:]))[:, 1].cpu()
            # take simple average
            consensus_est = (exp_seq_est +
                             exp_motif_est) / 2

            # calculate AUC for each
            seq_auc = metrics.roc_auc_score(labels, exp_seq_est)
            motif_auc = metrics.roc_auc_score(labels, exp_motif_est)
            consensus_auc = metrics.roc_auc_score(labels, consensus_est)

            # store model if perfomance is better than current best
            if seq_auc > prev_seq_auc:
                seq_state = sequence_model.state_dict()
                prev_seq_auc = seq_auc

            if motif_auc > prev_motif_auc:
                motif_state = motif_model.state_dict()
                prev_motif_auc = motif_auc

            # print out performance
            print("Test Set Results:")
            print(f'Sequence Model AUC: {seq_auc}')
            print(f'Motif Model AUC: {motif_auc}')
            print(f'Consensus Model AUC: {consensus_auc}')
            print("\n")

    # return to train mode for next iteration
    sequence_model.train()
    motif_model.train()


# look at ultimate performance
t_dat, t_labels = next(iter(test_loader))
# reset model states to best performance
sequence_model.load_state_dict(seq_state)
sequence_model.eval()
motif_model.load_state_dict(motif_state)
motif_model.eval()

# performance with seq only was best so use to determine performance
seq_est = torch.exp(sequence_model(t_dat.type(dtype)[:, :, :20]))[:, 1].cpu()
motif_est = torch.exp(motif_model(t_dat.type(dtype)[:, :, 22:]))[:, 1].cpu()

overall_est = (seq_est + motif_est)/2

preds_out = pd.DataFrame(list(zip(t_labels, overall_est)), columns=['class',
                                                                    'pred'])
preds_out.to_csv("./data/validation_data/output/epitope_test_" +
                 str(options.window_size)+"aa_preds.csv", index=False)

# call classes so that sensitivity and specificity can be calculated
seq_guess_class = []
for est in seq_est:
    if est >= .5:  # changing threshold alters se and sp, .5 = default
        seq_guess_class.append(1)
    else:
        seq_guess_class.append(0)

# print AUC
seq_auc = metrics.roc_auc_score(t_labels.detach().numpy(),
                                seq_est.detach().numpy())
print(seq_auc)

# Print classification report
seq_report = metrics.classification_report(t_labels.detach().numpy(),
                                           seq_guess_class)
print(seq_report)

# calculate and print se and sp
tn, fp, fn, tp = metrics.confusion_matrix(t_labels.detach().numpy(),
                                          seq_guess_class).ravel()
sensitivity = tp/(tp + fn)
specificity = tn/(tn+fp)

print("Sensitivity: ", sensitivity)
print("Specificity: ", specificity)

if options.human_only:
    torch.save(seq_state, options.out + "/human_epitope_sequence_mod.pt")
    torch.save(motif_state, options.out + "/human_epitope_motif_mod.pt")
else:
    torch.save(seq_state, options.out + "/all_mammal_epitope_sequence_mod.pt")
    torch.save(motif_state, options.out + "/all_mammal_epitope_motif_mod.pt")


# export saliencies
saliency = Saliency(motif_model)

pos_saliency_loader = torch.utils.data.DataLoader(pos_train,
                                                  batch_size=len(pos_train))
pos_saliency = next(iter(pos_saliency_loader))

print(len(pos_saliency))
print(pos_saliency.shape)

pos_grads = saliency.attribute(pos_saliency[:, :, 22:].type(dtype),
                               target=1,
                               abs=True)

neg_saliency_loader = torch.utils.data.DataLoader(neg_train,
                                                  batch_size=len(neg_train))

neg_saliency = next(iter(neg_saliency_loader))
neg_grads = saliency.attribute(neg_saliency[:, :, 22:].type(dtype),
                               target=0,
                               abs=True)

pos_mean = pos_grads.mean(axis=0)
neg_mean = neg_grads.mean(axis=0)

motif_df = pd.DataFrame(pos_mean, dtype=float)
motif_df = motif_df.append(pd.DataFrame(neg_mean, dtype=float))

motif_df['cleave_type'] = ['y']*len(pos_mean) + ['n']*len(neg_mean)


saliency = Saliency(sequence_model)
pos_grads = saliency.attribute(pos_saliency[:, :, :20].type(dtype),
                               target=1,
                               abs=True)

neg_grads = saliency.attribute(neg_saliency[:, :, :20].type(dtype),
                               target=0,
                               abs=True)

pos_mean = pos_grads.mean(axis=0)
neg_mean = neg_grads.mean(axis=0)

sequence_df = pd.DataFrame(pos_mean, dtype=float)
sequence_df = sequence_df.append(pd.DataFrame(neg_mean, dtype=float))

sequence_df['cleave_type'] = ['y']*len(pos_mean) + ['n']*len(neg_mean)


# motif_df.to_csv("./data/validation_data/output/plots/epitope_motif_saliency.csv", index=False)
# sequence_df.to_csv("./data/validation_data/output/plots/epitope_sequence_saliency.csv", index=False)
