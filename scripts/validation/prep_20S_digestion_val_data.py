#!/usr/bin/env python3
"""
prep_20S_digestion_val_data.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script loads 20S validation data and remaps column headers for consistency
"""

import pandas as pd
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--in-file",
                  help="pickled dictionary of cleavage windows")
parser.add_option("-o", "--out",
                  help="output directory where results will be exported")

(options, args) = parser.parse_args()

digestion_val = pd.read_csv(options.in_file)

new_digestion_cols = ['lit_reference', 'protein_name', 'origin_species',
                      'Proteasome', 'Subunit', 'full_seq_accession', 'end_pos',
                      'fragment', 'full_sequence', 'start_pos', 'entry_source']

digestion_val.columns = new_digestion_cols
digestion_val['exclusions'] = None

digestion_val.to_csv(options.out_file, index=False)
