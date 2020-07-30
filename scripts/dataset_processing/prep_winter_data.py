#!usr/bin/env python3
"""
prep_winder_data.py

For issues contact Ben Weeder (weeder@ohsu.edu)
"""

import pandas as pd
from extraction_functions import *
from optparse import OptionParser


# define command line parameters
parser = OptionParser()
parser.add_option("-i", "--in_file", dest="in_file",
                  help="CSV of data from winter et al.")
parser.add_option("-o", "--out_dir", dest="out_dir",
                  help="output directory where processed CSV is exported")

(options, args) = parser.parse_args()


dat = pd.read_csv(options.in_file, low_memory=False, header=0)
col_names = ['fragment', 'start_pos', 'end_pos', 'full_sequence',
             'exclusions', 'Proteasome']

out_df = pd.DataFrame(columns=col_names)

for i in range(len(dat)):
    row = dat.iloc[i]
    fragment_entries = parse_cleavage_logo(row)
    for fragment in fragment_entries:
        new_entry = pd.Series(fragment, index=col_names)
        out_df = out_df.append(new_entry, ignore_index=True)

out_df['Name'] = "NA"
out_df['Organism'] = "human"
out_df['Subunit'] = "20S"
out_df['entry_source'] = "cleavage_map"
out_df['DOI'] = "10.7554/eLife.27364"

out_df.to_csv(options.out_dir + "/winter_et_al_cleavage_fragments.csv",
              index=False)
