#!usr/bin/env python3
"""
pepsickle_pt_data_processing.py

For issues contact Ben Weeder (weeder@ohsu.edu)

pulls Cleavage information from pt. epitopes
"""

import pandas as pd
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--in-vivo",
                  help="input (epitope) file...")
parser.add_option("-b", "--in-vitro",
                  help="input file 2...")
parser.add_option("-o", "--out")

(options, args) = parser.parse_args()

dat = pd.read_csv(options.in_vivo, sep="\t")
dat_iv = pd.read_csv(options.in_vitro, sep="\t")

dat.columns = [x.strip() for x in dat.columns]
dat_iv.columns = [x.strip() for x in dat_iv.columns]
unique_ids = dat['protein_id'].unique()

term_cleave = []
term_prob = []
upstream_cleave = []
upstream_cleave_probs = []
max_upstream = []
min_upstream = []
response = []
internal_num = []

for p_id in unique_ids:
    subset = dat[dat['protein_id'] == p_id]
    subset_iv = dat_iv[dat_iv['protein_id'] == p_id]
    term_cleave.append(subset['cleaved'].iloc[-9].item())
    term_prob.append(subset['cleav_prob'].iloc[-9].item())
    upstream = True in list(subset_iv['cleaved'][8:16])
    upstream_cleave.append(upstream)
    # for upstream only
    tmp_cleave_probs = list(subset_iv['cleav_prob'][8:16])

    # for full window
    # tmp_cleave_probs = list(subset_iv['cleav_prob'][0:-9])
    # tmp_cleave_probs[-9] = subset['cleav_prob'].iloc[-9].item()

    upstream_cleave_probs.append(1 - np.prod(list(tmp_cleave_probs)))
    max_upstream.append(max(subset_iv['cleav_prob'][8:16]))
    min_upstream.append(min(subset_iv['cleav_prob'][8:16]))

    # num_cleaved = sum(list(subset_iv['cleaved']))/len(list(subset_iv['cleaved']))
    internal_cleave = 1 - np.prod(list(subset_iv['cleav_prob'][16:-9]))

    internal_num.append(internal_cleave)
    # max_internal = max(list(subset_iv['cleav_prob'][16:-9]))
    min_internal = min(list(subset_iv['cleav_prob'][16:-9]))
    # internal_num.append(min_internal)

    response.append(p_id.split("_")[-1])


out_df = pd.DataFrame(zip(term_cleave, term_prob, upstream_cleave,
                          max_upstream, internal_num, response, upstream_cleave_probs, min_upstream),
                      columns=['term_cleaved', "term_prob", 'upstream_cleaved',
                               'max_upstream_prob', 'num_internal', 'response',
                               'upstream_probs', 'min_upstream'])

both_true = out_df['term_cleaved'] & out_df['upstream_cleaved']
print(pd.crosstab(both_true, out_df['response']))

out_df.to_csv(options.out)
