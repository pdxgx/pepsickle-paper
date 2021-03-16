#!usr/bin/env python3
"""
epitope_index_check.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script verifies the location of each reported fragment within its source
and will attempt to repair indices for misreported locations.

-i, --in_file: file with reported fragments and indices to be verified
-o, --out_file: file with all verified fragments
"""

from sequence_featurization_tools import *
import re
import pandas as pd
from optparse import OptionParser

# define command line parameters
parser = OptionParser()
parser.add_option("-i", "--in_file", dest="in_file",
                  help="pandas data frame of fully merged fragment info")
parser.add_option("-o", "--out_file", dest="out_file",
                  help="output file where vetted CSV is written")

(options, args) = parser.parse_args()

df = pd.read_csv(options.in_file, low_memory=False)
df = df.where(pd.notnull(df), None)
df = df[df.full_sequence.notnull()]
df = df.reset_index(drop=True)
print("Total unchecked entries: ", len(df))


failed_indices = []
for i in range(len(df)):
    # pull out a single row
    row = df.iloc[i]
    # convert the epitope into a regular expression for seq. search
    search_exp = create_sequence_regex(row['fragment'])
    # if there's a starting and ending position...
    if row['start_pos'] and row['end_pos']:
        # get start & stop, these may need shifted depending on ref. indexing
        start = int(row['start_pos']) - 1  # if base 0
        end = int(row['end_pos'])
        # pull out expected match from source
        source_seq = str(row['full_sequence'][start:end])
        # if the sequences match, go to next row iteration
        if re.match(search_exp, source_seq):
            df.at[i, 'start_pos'] = start
            continue
        # if doesn't match, but only occurs once in origin sequence...
        elif len(re.findall(search_exp, str(row['full_sequence']))) == 1:
            # find location
            s = re.search(search_exp, str(row['full_sequence']))
            # update start position with start of match
            df.at[i, 'start_pos'] = s.span()[0]
            # update end position with start of match
            df.at[i, 'end_pos'] = s.span()[1]
        # if 0 or >1 matches, fail
        else:
            failed_indices.append(i)
    # either start or end is missing but not both
    elif row['start_pos'] or row['end_pos']:
        # if start position is present
        if row['start_pos']:
            start = int(row['start_pos']) - 1  # if 0 base
            end = int(row['start_pos']) + \
                  len(row['fragment'])
        # if end position is present
        if row['end_pos']:
            end = int(row['end_pos'])
            start = end - len(row['fragment']) - 1  # if 0 base
        # with newly defined positions, repeat code block from above
        source_seq = str(row['fragment'][start:end])
        # if the sequences match, go to next row iteration
        if re.match(search_exp, source_seq):
            df.at[i, 'start_pos'] = start
            df.at[i, 'end_pos'] = end
            continue
        # if doesn't match, but only occurs once in origin sequence...
        elif len(re.findall(search_exp, str(row['full_sequence']))) == 1:
            # find location
            s = re.search(search_exp, str(row['full_sequence']))
            # update start position with start of match
            df.at[i, 'start_pos'] = s.span()[0]
            # update end position with start of match
            df.at[i, 'end_pos'] = s.span()[1]
        # if 0 or >1 matches, fail
        else:
            failed_indices.append(i)
    # if no start and no end position
    else:
        # if only one occurrence...
        if len(re.findall(search_exp, str(row['full_sequence']))) == 1:
            # find location
            s = re.search(search_exp, str(row['full_sequence']))
            # update start position with start of match
            df.at[i, 'start_pos'] = s.span()[0]
            # update end position with start of match
            df.at[i, 'end_pos'] = s.span()[1]
        # if 0 or >1 matches, fail
        else:
            failed_indices.append(i)

print("Number of failed indices: ", len(failed_indices))
df.drop(index=failed_indices, inplace=True)
df = df.reset_index(drop=True)

# double check...
mismatch_indices = []
for i in range(len(df)):
    row = df.iloc[i]
    start = int(row['start_pos'])
    end = int(row['end_pos'])
    source_seq = str(row['full_sequence'][start:end])
    if str(row['fragment']) == source_seq:
        pass
    else:
        mismatch_indices.append(i)

# for now, drop few with errors:
print("Number of mismatch indices: ", len(mismatch_indices))
df.drop(index=mismatch_indices, inplace=True)
df = df.reset_index(drop=True)

print("Total entries left: ", len(df))

print("Data Summary: ")
print(df['Subunit'].value_counts())

df_20S = df[df['Subunit'] == "20S"]
if len(df_20S) > 0:
    print("20S Total Unique Fragments: ",
          len(df_20S['fragment'].unique()))
    print("20S Total Unique Source proteins: ",
          len(df_20S['full_sequence'].unique()))
    print("20S Fragment to Source Ratio: ",
          len(df_20S['fragment'].unique())/len(df_20S['full_sequence'].unique()))

df_26S = df[df['Subunit'] == "26S"]
if len(df_26S) > 0:
    print("26S Total Unique Source proteins: ",
          len(df_26S['full_sequence'].unique()))
    print("26S Fragment to Source Ratio: ",
          len(df_26S['fragment'].unique())/len(df_26S['full_sequence'].unique()))

df_epitope = df[df['Subunit'].isnull()]
if len(df_epitope) > 0:
    print("Epitope unique fragments: ",
          len(df_epitope['fragment'].unique()))
    print("Epitope Total Unique Source proteins: ",
          len(df_epitope['full_sequence'].unique()))
    print("Epitope Fragment to Source Ratio: ",
          len(df_epitope['fragment'].unique())/len(df_epitope['full_sequence'].unique()))

df.to_csv(options.out_file, index=False)
