#!/usr/bin/env python3
"""
prep_epitope_validation_data.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script extracts MHC class I epitope examples from [DOI] and maps them to
their source sequences. results are exported to a CSV and used downstream for
model validation_prep
"""

from optparse import OptionParser
import torch
import os
import pickle

parser = OptionParser()
parser.add_option("-i", "--in-dir",
                  help="directory of pytorch model weights")
parser.add_option("-o", "--out",
                  help="output directory where summary dict will be exported")

(options, args) = parser.parse_args()

model_dict = {}
for f in os.listdir(options.in_dir):
    if not f.startswith("."):
        if f.endswith(".pt"):
            path = options.in_dir + "/" + f
            mod_name, ext = os.path.splitext(f)
            mod = torch.load(path, map_location=torch.device('cpu'))
            model_dict[mod_name] = mod

out_file = options.out + "/trained_model_dict.pickle"
pickle.dump(model_dict, open(out_file, 'wb'))
