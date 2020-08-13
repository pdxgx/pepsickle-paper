#!usr/bin/env python3
"""

"""

import pandas as pd
from sequence_featurization_tools import *
from optparse import OptionParser

# define command line parameters
parser = OptionParser()
parser.add_option("-i", "--in-file",
                  help="CSV of data from levy et al.")
parser.add_option("-o", "--out",
                  help="output directory where processed CSV is exported")

(options, args) = parser.parse_args()

dat = pd.read_csv(options.in_file, low_memory=False)
unique_fragment_ids = dat['frag_tracker'].unique()

out_df = pd.DataFrame(columns=dat.columns)

for f_id in unique_fragment_ids:
    tmp_dat = dat[dat['frag_tracker'] == f_id]
    upstream_list = []
    downstream_list = []
    for _, entry in tmp_dat.iterrows():
        up_window = get_peptide_window(entry['full_sequence'],
                                       entry['start_pos'], entry['start_pos'])
        upstream_list.append(up_window)

        down_window = get_peptide_window(entry['full_sequence'],
                                         entry['end_pos'], entry['end_pos'])
        downstream_list.append(down_window)
    if (all(x == upstream_list[0] for x in upstream_list) and
            all(x == downstream_list[0] for x in downstream_list)):
        out_df = out_df.append(entry)
    else:
        continue

out_df.to_csv(options.out, index=False)
