#!usr/bin/env python3
"""
process_pepsickle_results.py

For issues contact Ben Weeder (weeder@ohsu.edu)

this script processes the output of pepsickle and returns summary metrics for
comparison with other cleavage tools.
"""

import pandas as pd
from sklearn import metrics
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-i", "--in-file", dest="in_file",
                  help="input file...")
parser.add_option("-o", "--out")

(options, args) = parser.parse_args()

dat = pd.read_csv(options.in_file, sep="\t")
dat.columns = [i.strip() for i in dat.columns]
unique_ids = list(dat['protein_id'].unique())

true_labels = []
prob_list = []
pred_list = []

for p_id in unique_ids:
    sub_dat = dat[dat['protein_id'] == p_id]
    mid_point = int(len(sub_dat)/2)

    prob_list.append(sub_dat['cleav_prob'].iloc[mid_point])
    pred_list.append(sub_dat['cleaved'].iloc[mid_point])

    if "pos" in str(p_id):
        true_labels.append(1)
    else:
        true_labels.append(0)

out_df = pd.DataFrame(list(zip(unique_ids, true_labels, prob_list, pred_list)),
                      columns=['ID', 'true_label', 'cleavage_prob',
                               'cleavage_pred'])
out_df.to_csv(options.out)
