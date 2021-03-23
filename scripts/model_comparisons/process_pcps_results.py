#!usr/bin/env python
"""
process_pcps_results.py

For issues contact Ben Weeder (weeder@ohsu.edu)

this script processes the output of PCPS and returns summary metrics for
comparison with other cleavage tools.
"""

from __future__ import print_function
from sklearn.metrics import (classification_report, roc_auc_score)
import argparse
import os
import sys
import pandas as pd


def get_predictions(site):
	''' Get cleavage model predictions

		site: list containing scores for the putative cleavage site;
			position 1 is for constitutive proteasome model, position 
			3 is for immunproteasome model

		Return value: constitutive model prediction, immuno model prediction
	'''
	if site[1] == 'YES' and site[3] == 'YES':
		return 1, 1
	elif site[1] == 'YES':
		return 1, 0
	elif site[3] == 'YES':
		return 0, 1
	else:
		return 0, 0


if __name__ == '__main__':

	# Set up command line parameters
	parser = argparse.ArgumentParser()
	parser.add_argument(
				'--input-file', '-i', type=str, required=True,
				help='path to input file with PCPS results'
	)
	parser.add_argument(
				'--peptide-size', '-p', type=int, required=True,
				help='size of peptides being assessed'
	)
	parser.add_argument(
				'--out', '-o', type=str, required=True,
				help='csv summary file.'
	)
	args = parser.parse_args()

	# Get path to input file
	input_file = os.path.abspath(args.input_file)
	output_file = os.path.abspath(args.out)
	# Store list of labels and predictions
	true_labels = []
	constitutive_predictions = []
	immuno_predictions = []
	constitutive_probs = []
	immuno_probs = []
	identifiers = []

	# Get peptide size + putative cleavage position (0-based)
	size = args.peptide_size
	position = int((size-1)/2)

	# Get predictions
	label = None
	scores = []
	with open(input_file) as f:
		for line in f:
			if '<br' in line:
				# Skip warning lines
				continue
			if line.startswith(' >') or line.startswith('>'):
				# New sequence
				if label is not None:
					assert len(scores) == size
					site = scores[position]
					# Get constitutive model prediction
					cons, immuno = get_predictions(site)
					constitutive_predictions.append(cons)
					immuno_predictions.append(immuno)
					cons_p = site[2]
					immuno_p = site[4]
					constitutive_probs.append(cons_p)
					immuno_probs.append(immuno_p)
					identifiers.append(identifier)
				# Get new label and scores list
				scores = []
				if line.startswith(' > N E G') or line.startswith('> N E G'):
					label = 0
					identifier = line.split(",")[0].replace(">", "").replace(" ", "")
				else:
					label = 1
					identifier = line.split(",")[0].replace(">", "").replace(" ", "")
				true_labels.append(label)
			elif not line.startswith('Aminoacid'):
				scores.append(line.strip().split(','))
	# Grab final predictions
	assert len(scores) == size
	site = scores[position]
	cons, immuno = get_predictions(site)

	constitutive_predictions.append(cons)
	immuno_predictions.append(immuno)

	constitutive_probs.append(cons_p)
	immuno_probs.append(immuno_p)

	identifiers.append(identifier)

	# Get performance evaluations
	print('Constitutive Model Results:')
	constitutive_report = classification_report(true_labels, constitutive_predictions)
	print(constitutive_report)
	constitutive_auc = roc_auc_score(true_labels, constitutive_predictions)
	print('AUC:', constitutive_auc)
	print('Immuno Model Results:')
	immuno_report = classification_report(true_labels, immuno_predictions)
	print(immuno_report)
	immuno_auc = roc_auc_score(true_labels, immuno_predictions)
	print('AUC:', immuno_auc)

	out_df = pd.DataFrame(list(zip(identifiers, true_labels,
							  constitutive_predictions,
							  constitutive_probs,
							  immuno_predictions,
							  immuno_probs)),
						  columns=["ID", "true_label", "constit_mod_pred",
								   "constit_prob", "immuno_prob_pred",
								   "immuno_prob"])
	out_df.to_csv(output_file)
