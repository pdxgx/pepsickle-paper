#!/usr/bin/env python
"""
proteasome_GBM_model.py

For issues contact Ben Weeder (weeder@ohsu.edu)

This script trains a gradient boosted model on the specified input features and
exports the best performing model (test AUC) as a model.joblib file
"""

from __future__ import print_function
from datetime import datetime
from joblib import dump
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.metrics import (
								classification_report, 
								roc_auc_score, roc_curve

)
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import random
import sys

# From the https://github.com/pdxgx/proteasome/ repository
# Located in data_processing/scripts/sequence_featurization_tools.py
_features = {
	'A': [
			1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,    6.0,   56.15265, -0.495, -2.4
		 ],
	'C': [
			0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   5.07,   69.61701,  0.081, -4.7
		 ],
	'D': [
			0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   2.77,   70.04515,  9.573, -4.5
		 ],
	'E': [
			0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   3.22,   86.35615,  3.173, -5.2
		 ],
	'F': [
			0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1,      0,   5.48,  119.722, -0.370, -4.9
		 ],
	'G': [
			0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   5.97,    37.80307,  0.386, -1.9
		 ],
	'H': [
			0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1,      0,   7.59,   97.94236,  2.029, -4.4
		 ],
	'I': [
			0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   6.02,  103.6644, -0.528, -6.6
		 ],
	'K': [
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   9.74,   102.7783,  2.101, -7.5
		 ],
	'L': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   5.98,  102.7545, -0.342, -6.3
		 ],
	'M': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   5.74,  103.928, -0.324, -6.1
		 ],
	'N': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   5.41,   76.56687,  2.354, -4.7
		 ],
	'P': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
			0,      0,    6.3,   71.24858, -0.322, -0.8
		 ],
	'Q': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			0,      0,   5.65,   88.62562,  2.176, -5.5
		 ],
	'R': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
			0,      0,  10.76,  110.5867,  4.383, -6.9
		 ],
	'S': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
			0,      1,   5.68,   55.89516,  0.936, -4.6
		 ],
	'T': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
			0,      1,    5.6,   72.0909,  0.853, -5.1
		 ],
	'V': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
			0,      0,   5.96,   86.28358, -0.308, -4.6
		 ],
	'W': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
			1,      0,   5.89,  137.5186,  -0.27, -4.8
		 ],
	'Y': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
			1,      1,   5.66,  121.5862,  1.677, -5.4
		 ],
	'*': [
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0, 7.5,      0.0, 1.689157, 0.0
		 ],
	'B': [
			0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   4.09,   73.30601,  5.964, -4.6
		 ],
	'Z': [
			0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			0,      0,   4.44,   87.49089,  2.675, -5.35
		 ],
	'J': [
			0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,    6.0,  103.2094,  -0.426, -6.45
		 ],
	'U': [
			0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0,      0,   5.07,   69.61701,  0.081, -4.7
		 ],
	'X': [
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			0.5, 0.5,    6.008095,    88.55829,    0.6195, -4.845
		 ]
}


def time_print(string):
	""" Prints string with date/time

		string: string to print

		No return value.
	"""
	print(datetime.now(), string)


def split_range(string, descriptor):
	""" Splits range string into list of values

		string: string to split
		descriptor: description of option for error message

		Return value: list of values
	"""
	if ',' in string:
		# Get list of values
		values = [int(x) for x in string.split(',')]
		try:
			assert len(values) == 3
			return [
						int(x) for x in np.linspace(
														start = int(values[0]), 
														stop = int(values[1]), 
														num = int(values[2])
										)
					]
		except AssertionError:
			raise ValueError(
								''.join(
											[
												'Range for ', descriptor,
												' may be given in threes ',
												'(min,max,count) or as a ',
												'single number'
											]
										)
			)
	else:
		# Create list of one value
		return [int(string)]

# From the https://github.com/pdxgx/proteasome/ repository
# Located in data_processing/scripts/sequence_featurization_tools.py
def featurize_sequence(seq, normalize=False):
	"""
	takes an input aa sequence of any length and returns a 2D numpy array
	of feature values with rows as positions and columns as feature values
	:param seq: a string of amino acid symbols of any length
	:return feature_matrix:
	"""
	feature_matrix = np.array([_features[aa] for aa in seq], dtype=float)
	return feature_matrix

def extended_featurization(seq, annotations, normalize=False):
	""" Creates feature array and adds binary constitutive/immunoproteasome info

		seq: sequence to featurize
		annotation: set of annotations for seq, with proteasome type in position
					5 of each annotation

		Return value: list of features for sequence
	"""
	featurized = featurize_sequence(seq).tolist()
	proteasome = [x[1] for x in annotations]
	if 'C' in proteasome or 'M' in proteasome:
		featurized.append(1)
	else:
		featurized.append(0)
	if 'I' in proteasome or 'M' in proteasome:
		featurized.append(1)
	else:
		featurized.append(0)
	return featurized


# Adapted from https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
def flatten(l):
	"""
	Flattens list of lists into one list

	list: list of lists
	
	Return value: flattened list
	"""
	return [item for sublist in l for item in sublist]


def evaluate(classifier, features, labels, plot_file=None,
			 csv_file=None):
	""" 
	Evaluates classifier using accuracy, precision, and recall scores;
	generates an ROC curve to display or save.

	classifier: classifier model
	features: feature matrix used to make predictions for each sample
	labels: true labels for each sample

	Return values: accuracy, precision, recall; prints or saves ROC plot
	"""
	# Predict labels
	predictions = classifier.predict(features)
	# Get evaluation scores
	report = classification_report(labels, predictions)
	print(report)
	# Generate data for ROC curve
	probs = classifier.predict_proba(features)
	fpr, tpr, thresholds = roc_curve(labels, probs[:, 1])
	auc = roc_auc_score(labels, probs[:, 1])
	print('AUC:', auc)
	# Plot ROC curve
	plt.plot(
				fpr, tpr, color='blueviolet', lw=2, 
				label='ROC curve (AUC = %0.2f)' % auc
	)
	plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
	plt.xlim([-0.05, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate (FPR)')
	plt.ylabel('True Positive Rate (TPR)')
	plt.legend(loc="lower right")
	# Save to plot to file or display
	if plot_file is not None:
		plt.savefig(plot_file)
		plt.close()
	else:
		plt.show()
		plt.close()

	if csv_file:
		out_df = pd.DataFrame(zip(labels, probs[:, 1], predictions),
							  columns=['true_label', 'cleave_p', 'cleave_pred'])
		out_df.to_csv(csv_file)


if __name__ == '__main__':

	# Set up command line parameters
	parser = argparse.ArgumentParser()
	parser.add_argument(
				'--input-file', '-i', type=str, 
				default='/Users/wooma/proteasome_window_size_6.pickle',
				help='path to pickled dictionary with input data'
	)
	parser.add_argument(
				'--output-directory', '-o', type=str, 
				default=None,
				help='path to directory to store output files'
	)
	parser.add_argument(
				'--estimators', '-e', type=str, default="100,1000,10",
				help='range for number of estimators to use: min,max,count'
					 ' OR number'
	)
	parser.add_argument(
				'--tuning-combinations', '-t', type=int, default=100,
				help='number of hyperparameter combinations to tune with'
	)
	parser.add_argument(
				'--cross-validation', '-k', type=int, default=10, # was 4
				help='number of folds for cross validation'
	)
	parser.add_argument(
				'--training-proportion', '-p', type=float, default=0.8,
				help='proportion of data to use for training'
	)
	parser.add_argument(
				'--model-type', '-m', type=str, default='epitope',
				help='model type to generate (epitope or proteasome)'
	)
	parser.add_argument(
				'--features', '-f', type=str, default='identity',
				help='features to use (identity, chemistry, or both)'
	)
	parser.add_argument(
				'--normalize', '-n', required=False, action='store_true',
				help='whether to normalize chemical features'
	)
	args = parser.parse_args()

	time_print('Loading and preprocessing input data...')

	# Set random seed
	# seed_num = 1991  # was 427
	random.seed(427)  # was 427
	rng = np.random.RandomState(1991)  # was 423

	# Import pickled dictionary
	with open(os.path.abspath(args.input_file), 'rb') as p:
		data_dict = pickle.load(p)

	if args.normalize:
		# Normalize values in _features dictionary
		min_max_scaler = MinMaxScaler()
		df = pd.DataFrame.from_dict(_features, orient='index')
		minmax_df = pd.DataFrame(min_max_scaler.fit_transform(df), index=df.index)
		for index, row in minmax_df.iterrows():
			for i in range(22, len(row)):
				_features[index][i] = minmax_df.loc[index, i]

	# Shuffle labeled data, converting peptide windows to feature arrays
	if args.model_type == 'epitope':
		labeled_data = [(featurize_sequence(x).tolist(), 1) for x in data_dict['epitope']['positives']]
		labeled_data.extend([(featurize_sequence(x).tolist(), 0) for x in data_dict['epitope']['negatives']])
	else:
		labeled_data = [(extended_featurization(x, data_dict['proteasome']['positives'][x]), 1) for x in data_dict['proteasome']['positives']]
		labeled_data.extend([(extended_featurization(x, data_dict['proteasome']['negatives'][x]), 0) for x in data_dict['proteasome']['negatives']])
	random.shuffle(labeled_data)

	c_p = []
	c_n = []
	i_p = []
	i_n = []
	for entry in labeled_data:
		if entry[1] == 1:
			if entry[0][-1] == 1:
				i_p.append(1)
			if entry[0][-2] == 1:
				c_p.append(1)
		else:
			if entry[0][-1] == 1:
				i_n.append(1)
			if entry[0][-2] == 1:
				c_n.append(1)

	print("Total Constit pos: ", sum(c_p))
	print("Total Constit neg: ", sum(c_n))
	print("Total immuno pos: ", sum(i_p))
	print("Total immuno neg: ", sum(i_n))
	print("Total Unique: ", len(labeled_data))

	# Get hyperparameters to tune with
	n_estimators = split_range(args.estimators, '-e/--estimators')
	max_features = ['auto', 'sqrt', 'log2']
	# random_grid = {'n_estimators': n_estimators, 'max_features': max_features, 'class_weight': ['balanced']}
	random_grid = {'n_estimators': n_estimators, 'max_features': max_features}

	# Select feature indices:
	if args.features == 'identity':
		start = 0
		stop = 20
	elif args.features == 'chemistry':
		start = 22
		stop = 26
	else:
		start = 0
		stop = 26

	# Split training, testing, and validation data
	if args.model_type == 'epitope':
		X = [flatten([y[start:stop] for y in x[0]]) for x in labeled_data]
	else:
		X = [flatten([y[start:stop] for y in x[0] if y != x[0][-1] and y != x[0][-2]]) + list(x[0][-2:]) for x in labeled_data]
	y = [x[1] for x in labeled_data]
	X_train, X_test, y_train, y_test = train_test_split(X, y,
														train_size=args.training_proportion,
														random_state=rng)

	sample_weight_vals = []
	pos_prop = sum(y_train)/len(y_train)
	for entry in y_train:
		if entry == 1:
			sample_weight_vals.append(1 - pos_prop)
		else:
			sample_weight_vals.append(pos_prop)

	time_print('Training model and tuning hyperparameters...')

	# Train classifier
	# classifier = RandomForestClassifier()
	classifier = GradientBoostingClassifier(random_state=rng)
	classifier_random = RandomizedSearchCV(
									estimator = classifier, 
									param_distributions = random_grid, 
									n_iter = args.tuning_combinations,
									cv = args.cross_validation, 
									n_jobs=-1,
									verbose=2,
									random_state=rng
	)

	# Fit model
	#classifier_random.fit(X_train, y_train)
	classifier_random.fit(np.array(X_train), np.array(y_train), sample_weight=sample_weight_vals)

	# Display best parameters from tuning
	print('Best parameters:')
	print(classifier_random.best_params_)
	
	time_print('Evaluating performance of classifier...')

	best_classifier = classifier_random.best_estimator_
	best_estimator_count = str(classifier_random.best_params_['n_estimators'])
	best_estimator_feats = str(classifier_random.best_params_['max_features'])

	feat_weights = best_classifier.feature_importances_
	max_weights = [(i, j) for i, j in enumerate(feat_weights) if j == max(feat_weights)]
	print('Most important feature(s), weight(s):', max_weights)

	# Get files for ROC plots/feature weights if applicable
	if args.output_directory is not None:
		train_roc_file = os.path.join(
										os.path.abspath(args.output_directory),
										'_'.join(['train_roc.rf', best_estimator_count, best_estimator_feats, '.png'])
		)
		test_roc_file = os.path.join(
										os.path.abspath(args.output_directory),
										'_'.join(['test_roc.rf', best_estimator_count, best_estimator_feats, '.png'])
		)
		weights_file = os.path.join(
										os.path.abspath(args.output_directory),
										'_'.join(['feature_weights.rf', best_estimator_count, best_estimator_feats, '.tsv'])
		)
		with open(weights_file, 'w') as f:
			print('\t'.join(['Feature', 'Weight']), file=f)
			for i in range(len(feat_weights)):
				print('\t'.join([str(i), str(feat_weights[i])]), file=f)
		train_file = os.path.join(
										os.path.abspath(args.output_directory),
										'_'.join(['training_set.tsv'])
		)
		train_frame = pd.DataFrame(X_train)
		train_frame['Y'] = y_train
		train_frame.to_csv(train_file, sep='\t', index=False)
		train_file = os.path.join(
										os.path.abspath(args.output_directory),
										'_'.join(['training_set.tsv'])
		)
		test_file = os.path.join(
										os.path.abspath(args.output_directory),
										'_'.join(['testing_set.tsv'])
		)
		test_csv = os.path.join(
										os.path.abspath(args.output_directory),
										'_'.join(['final_test_preds.csv']))

		test_frame = pd.DataFrame(X_test)
		test_frame['Y'] = y_test
		test_frame.to_csv(test_file, sep='\t', index=False)
	else:
		train_roc_file, val_roc_file, test_roc_file = None, None, None

	# calibrate
	# classifier_calibrated = CalibratedClassifierCV(best_classifier, cv="prefit")
	# classifier_calibrated.fit(X_test, y_test)
	classifier_calibrated = best_classifier
	# Evaluate on training, validation, and testing sets
	time_print('Evaluating performance on training set...')
	evaluate(classifier_calibrated, X_train, y_train, train_roc_file)
	time_print('Evaluating performance on testing set...')
	evaluate(classifier_calibrated, X_test, y_test, test_roc_file, csv_file=test_csv)

	dump(classifier_calibrated, os.path.join(args.output_directory, 'model.joblib'))

	# Plot feature importances
	feature_importance = permutation_importance(classifier_calibrated, X_test,
												y_test, n_repeats=30,
												random_state=0)
	# print(feature_importance.importances.shape)
	importance_df = pd.DataFrame(zip(list(feature_importance.importances_mean),
								  list(feature_importance.importances_mean)),
								 columns=['mean_importance', 'importance_stdev'])
	importance_df.to_csv(args.output_directory +
						 "/final_feature_importances.tsv", sep='\t',
						 index=False)

	print(datetime.now(), 'Done...')

	sys.exit()


