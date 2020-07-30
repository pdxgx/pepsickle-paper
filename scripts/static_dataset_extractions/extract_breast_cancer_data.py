#!/usr/bin/env Python3
"""extract_breast_cancer_data.py

Ellysia Li (lie@ohsu.edu)

Python 3.7

This script converts epitope information from a xlsx file into a csv file.
The information was derived from a study on breast cancer cell lines
(doi: 10.1016/j.jprot.2018.01.004)

This script requires that `pandas` and `numpy` be installed within the Python
environment you are running this script in.

Inputs:
    The location of the xlsx file of interest
    The location of a csv file containing two columns:
        ["Input"]: A column containing all entry names of UniProt proteins
                   which were derived from the constructed dataframe's
                   ["Protein Link"]
        ["Entry"]: The column of corresponding accession number (found using
                   UniProt's Retrieve/ID mapping function)

Outputs:
    The corresponding csv file containing the epitope's description and its
    protein's UniProt IRI
"""
import pandas as pd
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--in_dir", dest="in_dir",
                  help="Location of breast cancer epitope files")

parser.add_option("-o", "--out_dir", dest="out_dir",
                  help="output directory where CSV is exported")

(options, args) = parser.parse_args()
df = pd.DataFrame()

""" Data Extraction """
# the range(2, 24) corresponds to the different excel sheets which contain
# epitope information
base_file = pd.ExcelFile(options.in_dir + "/breast_cancer.xlsx")
for i in range(2, 24):
    tmp_df = pd.read_excel(options.in_dir + "/breast_cancer.xlsx", sheet_name=i,
                           header=2)[["Sequence", "Protein Link", "Unique"]]
    tmp_df['cell_line'] = base_file.sheet_names[i]
    df = df.append(tmp_df)

key = pd.read_csv(options.in_dir + "/breast_cancer_key.csv")
dict_key = {k: v for k, v in zip(key["Input"].to_list(),
                                 key["Entry"].tolist())}
df["Protein_ref"] = df["Protein Link"].map(dict_key)
df['Ref_type'] = "Uniprot"
df.dropna(subset=["Sequence"], inplace=True)
print("N entries: ", df.shape[0])

df.drop_duplicates(inplace=True)
print("N unique entries: ", df.shape[0])

df["Sequence"] = df["Sequence"].apply(lambda x: x.split(".")[1])

df.drop(columns=["Protein Link", "Unique"], inplace=True)
df.rename(columns={"Sequence": "fragment"}, inplace=True)

df.to_csv(options.out_dir + "/breast_cancer_epitopes.csv", index=False)
