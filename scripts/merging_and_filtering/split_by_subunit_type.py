#!usr/bin/env python3
"""
split_by_subunit_type.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script takes in processed cleavage fragments from degestion studies and devides them based on
cleavage subunit-complex type.
"""

import pandas as pd
from optparse import OptionParser


# define command line parameters
parser = OptionParser()
parser.add_option("-i", "--in-file",
                  help="CSV of data.")
parser.add_option("-o", "--out",
                  help="output directory where processed CSV is exported")
parser.add_option("--human-only", dest="human_only", action="store_true",
                  help='names exports with human only flag')

(options, args) = parser.parse_args()

# load data
dat = pd.read_csv(options.in_file, low_memory=False)

dat_20S = dat[dat['Subunit'] == '20S']
dat_26S = dat[dat['Subunit'] == '26S']
dat_epitope = dat[dat['Subunit'].isna()]

print("Total 20S entries", dat_20S.shape[0])
print("Total 26S entries", dat_26S.shape[0])
print("Total epitope entries", dat_epitope.shape[0])

if options.human_only:
    dat_20S.to_csv(options.out + "/merged_20S_fragments_human.csv",
                   index=False)
    dat_26S.to_csv(options.out + "/merged_26S_fragments_human.csv",
                   index=False)
    dat_epitope.to_csv(options.out + "/merged_epitope_fragments_human.csv",
                       index=False)
else:
    dat_20S.to_csv(options.out + "/merged_20S_fragments.csv", index=False)
    dat_26S.to_csv(options.out + "/merged_26S_fragments.csv", index=False)
    dat_epitope.to_csv(options.out + "/merged_epitope_fragments.csv",
                       index=False)
