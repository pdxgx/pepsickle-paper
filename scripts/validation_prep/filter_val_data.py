#!usr/bin/env python3
"""
filter_val_data.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script loads training and validation data to filter out examples preiously
seen during training and testing
"""

import os
import math
import pickle
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--in-dir",
                  help="directory with pickled dictionaries of validation_prep "
                       "cleavage windows")
parser.add_option("-t", "--training-data-dir",
                  help="directory with pickled dictionaries of training data")
parser.add_option("-w", "--window-size",
                  help="window size used for training and validation data")
parser.add_option("-o", "--out",
                  help="output directory where results will be exported")
parser.add_option("--internal-filter-size", default=None,
                  help="verifies center window of size N is unique for all "
                       "entries. If no value is specified, returns unique full"
                       "sized windows.")
parser.add_option("--human-only", dest="human_only", action="store_true",
                  default=False,
                  help="restricts to only human epitope exmaples in training data")

(options, args) = parser.parse_args()

if options.internal_filter_size:
    internal_size = int(options.internal_filter_size)
    buffer = (internal_size-1)/2
    center_point = math.ceil(int(options.window_size)/2)

in_file_list = os.listdir(options.in_dir)
if "20S_val_windows_" + str(options.window_size) + "aa_paired.pickle" in in_file_list:
    digestion_val_dict = pickle.load(
        open(options.in_dir + "/20S_val_windows_" + str(options.window_size) + "aa_paired.pickle", "rb"))
if "epitope_val_windows_" + str(options.window_size) + "aa_paired.pickle" in in_file_list:
    epitope_val_dict = pickle.load(
        open(options.in_dir + "/epitope_val_windows_" + str(options.window_size) + "aa_paired.pickle", "rb"))


training_files = os.listdir(options.training_data_dir)
if options.human_only:
    digestion_training_dict = pickle.load(
        open(options.training_data_dir +
             "/human_20S_windows_" + str(options.window_size) +
             "aa.pickle", "rb"))
    epitope_training_dict = pickle.load(
        open(options.training_data_dir +
             "/human_epitope_windows_" + str(options.window_size) +
             "aa.pickle", "rb"))

else:
    epitope_training_dict = pickle.load(
        open(options.training_data_dir +
             "/all_mammal_epitope_windows_" + str(options.window_size) +
             "aa.pickle", "rb"))
    digestion_training_dict = pickle.load(
        open(options.training_data_dir +
             "/all_mammal_20S_windows_" + str(options.window_size) + "aa.pickle", "rb"))


epitope_positive_windows = dict()
epitope_prev_neg_windows = dict()
for key in epitope_val_dict['epitope']['positives'].keys():
    if key not in epitope_training_dict['epitope']['positives'].keys():
        if key not in epitope_training_dict['epitope']['negatives'].keys():
            epitope_positive_windows[key] = epitope_val_dict['epitope']['positives'][key].copy()
        else:
            epitope_prev_neg_windows[key] = epitope_val_dict['epitope']['positives'][key].copy()


epitope_negative_windows = dict()
for key in epitope_val_dict['epitope']['negatives'].keys():
    if key not in epitope_training_dict['epitope']['positives'].keys():
        if key not in epitope_training_dict['epitope']['negatives'].keys():
            if key not in epitope_training_dict['epitope']['unknowns'].keys():
                epitope_negative_windows[key] = epitope_val_dict['epitope']['negatives'][key].copy()


if options.internal_filter_size:
    # track encountered internal windows
    tmp_epitope_negative_windows = epitope_negative_windows.copy()
    tmp_epitope_positive_windows = epitope_positive_windows.copy()
    prev_internal = []
    # perform positives first... this will filter negatives w/ matching internals
    for key in epitope_positive_windows.keys():
        internal_window = key[int(center_point-buffer):int(center_point+buffer)+1]
        if internal_window in prev_internal:
            tmp_epitope_positive_windows.pop(key)
        else:
            prev_internal.append(internal_window)

    for key in epitope_negative_windows.keys():
        internal_window = key[int(center_point-buffer):int(center_point+buffer)+1]
        if internal_window in prev_internal:
            tmp_epitope_negative_windows.pop(key)
        else:
            prev_internal.append(internal_window)

    print("Removed " +
          str((len(epitope_negative_windows)-len(tmp_epitope_negative_windows) +
               (len(epitope_positive_windows)-len(tmp_epitope_positive_windows)))) +
          " Epitope Windows")

    epitope_negative_windows = tmp_epitope_negative_windows.copy()
    epitope_positive_windows = tmp_epitope_positive_windows.copy()

epitope_windows_filtered = dict()
epitope_windows_filtered['positives'] = epitope_positive_windows
epitope_windows_filtered['negatives'] = epitope_negative_windows


digestion_constit_positive_windows = dict()
digestion_immuno_positive_windows = dict()
for key in digestion_val_dict['proteasome']['positives'].keys():
    if key not in digestion_training_dict['proteasome']['positives'].keys():
        if "C" in list(digestion_val_dict['proteasome']['positives'][key])[0]:
            digestion_constit_positive_windows[key] = digestion_val_dict['proteasome']['positives'][key].copy()

        if "I" in list(digestion_val_dict['proteasome']['positives'][key])[0]:
            digestion_immuno_positive_windows[key] = digestion_val_dict['proteasome']['positives'][key].copy()


digestion_constit_negative_windows = dict()
digestion_immuno_negative_windows = dict()
for key in digestion_val_dict['proteasome']['negatives'].keys():
    if key not in digestion_training_dict['proteasome']['positives'].keys():
        if key not in digestion_training_dict['proteasome']['negatives'].keys():
            if "C" in list(digestion_val_dict['proteasome']['negatives'][key])[0]:
                digestion_constit_negative_windows[key] = digestion_val_dict['proteasome']['negatives'][key].copy()
            if "I" in list(digestion_val_dict['proteasome']['negatives'][key])[0]:
                digestion_immuno_negative_windows[key] = digestion_val_dict['proteasome']['negatives'][key].copy()


if options.internal_filter_size:
    # track encountered internal windows
    tmp_digestion_constit_negative_windows = digestion_constit_negative_windows.copy()
    tmp_digestion_constit_positive_windows = digestion_constit_positive_windows.copy()
    prev_internal = []
    # perform positives first... this will filter negatives w/ matching internals
    for key in digestion_constit_positive_windows.keys():
        internal_window = key[int(center_point-buffer):int(center_point+buffer)+1]
        if internal_window in prev_internal:
            tmp_digestion_constit_positive_windows.pop(key)
        else:
            prev_internal.append(internal_window)

    for key in digestion_constit_negative_windows.keys():
        internal_window = key[int(center_point-buffer):int(center_point+buffer)+1]
        if internal_window in prev_internal:
            tmp_digestion_constit_negative_windows.pop(key)
        else:
            prev_internal.append(internal_window)

    print("Removed " +
          str((len(digestion_constit_negative_windows)-len(tmp_digestion_constit_negative_windows) +
               (len(digestion_constit_positive_windows)-len(tmp_digestion_constit_positive_windows)))) +
          " Constitutive Windows")

    digestion_constit_negative_windows = tmp_digestion_constit_negative_windows.copy()
    digestion_constit_positive_windows = tmp_digestion_constit_positive_windows.copy()

if options.internal_filter_size:
    # track encountered internal windows
    tmp_digestion_immuno_negative_windows = digestion_immuno_negative_windows.copy()
    tmp_digestion_immuno_positive_windows = digestion_immuno_positive_windows.copy()
    prev_internal = []

    # perform positives first... this will filter negatives w/ matching internals
    for key in digestion_immuno_positive_windows.keys():
        internal_window = key[int(center_point-buffer):int(center_point+buffer)+1]
        if internal_window in prev_internal:
            tmp_digestion_immuno_positive_windows.pop(key)
        else:
            prev_internal.append(internal_window)

    for key in digestion_immuno_negative_windows.keys():
        internal_window = key[int(center_point-buffer):int(center_point+buffer)+1]
        if internal_window in prev_internal:
            tmp_digestion_immuno_negative_windows.pop(key)
        else:
            prev_internal.append(internal_window)

    print("Removed " +
          str((len(digestion_immuno_negative_windows)-len(tmp_digestion_immuno_negative_windows) +
               (len(digestion_immuno_positive_windows)-len(tmp_digestion_immuno_positive_windows)))) +
          " Immuno Windows")

    digestion_immuno_negative_windows = tmp_digestion_immuno_negative_windows.copy()
    digestion_immuno_positive_windows = tmp_digestion_immuno_positive_windows.copy()


digestion_constit_windows_filtered = dict()
digestion_constit_windows_filtered['positives'] = digestion_constit_positive_windows
digestion_constit_windows_filtered['negatives'] = digestion_constit_negative_windows

digestion_immuno_windows_filtered = dict()
digestion_immuno_windows_filtered['positives'] = digestion_immuno_positive_windows
digestion_immuno_windows_filtered['negatives'] = digestion_immuno_negative_windows


# print summary stats
print("Epitope positives: ", len(epitope_positive_windows))
print("Epitope negatives: ", len(epitope_negative_windows))

print("Digestion constitutive positives: ", len(digestion_constit_positive_windows))
print("Digestion constitutive negatives: ", len(digestion_constit_negative_windows))

print("Digestion immuno positives: ", len(digestion_immuno_positive_windows))
print("Digestion immuno negatives: ", len(digestion_immuno_negative_windows))


# export dictionaries of validation_prep windows
pickle.dump(epitope_windows_filtered,
            open(options.out + "/window_dictionaries" +
                 "/epitope_val_filtered_" + str(options.window_size) + "aa.pickle", "wb"))

pickle.dump(digestion_constit_windows_filtered,
            open(options.out + "/window_dictionaries" +
                 "/20S_digestion_constitutive_validation_filtered_" + str(options.window_size) + "aa.pickle",
                 "wb"))

pickle.dump(digestion_immuno_windows_filtered,
            open(options.out + "/window_dictionaries" +
                 "/20S_digestion_immuno_validation_filtered_" + str(options.window_size) + "aa.pickle", "wb"))


# export as fasta files
val_handle = options.out + "/window_fasta_files"
epitope_val_fasta = open(val_handle + "/epitope_val_data_" + str(options.window_size) + "aa.fasta", "w")

for i in range(len(epitope_positive_windows)):
    prot_name = ">pos_" + str(i)
    epitope_val_fasta.write(prot_name)
    epitope_val_fasta.write("\n")
    epitope_val_fasta.write(list(epitope_positive_windows.keys())[i])
    epitope_val_fasta.write("\n")

for i in range(len(epitope_negative_windows)):
    prot_name = ">neg_" + str(i)
    epitope_val_fasta.write(prot_name)
    epitope_val_fasta.write("\n")
    epitope_val_fasta.write(list(epitope_negative_windows.keys())[i])
    epitope_val_fasta.write("\n")

epitope_val_fasta.close()


digestion_constit_val_fasta = open(val_handle +
                                   "/20S_digestion_constitutive_validation_data_" + str(options.window_size) + "aa.fasta", "w")

for i in range(len(digestion_constit_positive_windows)):
    prot_name = ">pos_" + str(i)
    digestion_constit_val_fasta.write(prot_name)
    digestion_constit_val_fasta.write("\n")
    digestion_constit_val_fasta.write(list(digestion_constit_positive_windows.keys())[i])
    digestion_constit_val_fasta.write("\n")

for i in range(len(digestion_constit_negative_windows)):
    prot_name = ">neg_" + str(i)
    digestion_constit_val_fasta.write(prot_name)
    digestion_constit_val_fasta.write("\n")
    digestion_constit_val_fasta.write(list(digestion_constit_negative_windows.keys())[i])
    digestion_constit_val_fasta.write("\n")

digestion_constit_val_fasta.close()


digestion_immuno_val_fasta = open(val_handle +
                                  "/20S_digestion_immuno_validation_data_" + str(options.window_size) + "aa.fasta",
                                  "w")

for i in range(len(digestion_immuno_positive_windows)):
    prot_name = ">pos_" + str(i)
    digestion_immuno_val_fasta.write(prot_name)
    digestion_immuno_val_fasta.write("\n")
    digestion_immuno_val_fasta.write(list(digestion_immuno_positive_windows.keys())[i])
    digestion_immuno_val_fasta.write("\n")

for i in range(len(digestion_immuno_negative_windows)):
    prot_name = ">neg_" + str(i)
    digestion_immuno_val_fasta.write(prot_name)
    digestion_immuno_val_fasta.write("\n")
    digestion_immuno_val_fasta.write(list(digestion_immuno_negative_windows.keys())[i])
    digestion_immuno_val_fasta.write("\n")

digestion_immuno_val_fasta.close()


