 #!/usr/bin/env python3
"""
prep_epitope_validation_data.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script extracts MHC class I epitope examples from
[https://doi.org/10.1038/ncomms13404] and maps them to their source sequences.
results are exported to a CSV and used downstream for model validation_prep
"""

import pandas as pd
from extraction_functions import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--in-file",
                  help="pickled dictionary of cleavage windows")
parser.add_option("-o", "--out",
                  help="output directory where results will be exported")

(options, args) = parser.parse_args()

# load data
dat = pd.read_csv(options.in_file, low_memory=False, skiprows=0, header=1)

# pull out class I columns
class_I_columns = []
for i in range(len(dat.columns)):
    col = dat.columns[i]
    if "_HLA-I " in col:
        class_I_columns.append(i)

# subset fragments presented on HLA class-I
total_class_I_intensity = dat[dat.columns[class_I_columns]].sum(axis=1)
dat = dat[total_class_I_intensity > 0]
# subset to relevant columns
cols_to_keep = ['Start position', 'End position', 'Sequence', 'Proteins']
dat = dat[cols_to_keep]
# drop incomplete examples
dat = dat.dropna()

# flag unambiguous mappings
unique_maps = []
for p_id in dat['Proteins']:
    id_count = p_id.count(";")
    if id_count == 0:
        unique_maps.append(True)
    else:
        unique_maps.append(False)

# subset to only those with unambiguous ID map
unambiguous_dat = dat[unique_maps]
unambiguous_dat.reset_index(inplace=True)

# get source sequences
source_sequences = {}
failed_entries = []
unique_prots = unambiguous_dat["Proteins"].unique()
for p in range(len(unique_prots)):
    p_id = unique_prots[p]
    try:
        entry = retrieve_UniProt_seq(p_id)
        source_sequences[p_id] = entry
    except:
        failed_entries.append(p_id)

    if p % 100 == 0:
        print(round(p/len(unique_prots), 3))

# map source sequences to examples
unambiguous_dat['full_sequence'] = None
for e in range(len(unambiguous_dat)):
    # pull protein ID
    prot_id = str(unambiguous_dat.loc[e, "Proteins"])
    # if sequence was found for given ID, store full sequence
    if prot_id in source_sequences.keys():
        unambiguous_dat.loc[e, 'full_sequence'] = source_sequences[prot_id]

# drop missing
unambiguous_dat = unambiguous_dat.dropna()
# rename columns
new_cols = ['index', 'start_pos',  'end_pos', 'fragment', 'full_seq_accession',
            'full_sequence']
unambiguous_dat.columns = new_cols
unambiguous_dat['entry_source'] = "validation_prep"

# export
unambiguous_dat.to_csv(options.out + "/validation_epitopes_w_source.csv",
                       index=False)
