# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022
## This script derives target binding profile of query compounds using the predictive MACCS fingerprints identified by TargetTox (feature selection pipeline)


## Module
import joblib
import numpy as np
import pandas as pd


## This function derives target binding profile of query compounds using the predictive MACCS fingerprints identified by TargetTox 
def derive_target_profile(structure_df):
	## 0. Input arguments:
		# structure_df: data frame contains MACSS fingprints of query compounds  
	
	## 1. Predict the target binding of query compounds using classifier built upon predictive fingerprints of each target 
	# read in TargetTox-identified MACCS fingerprints predictive of target binding   
	feature_file = 'data/binding/fingerprint_maccs_select_features_mc_0.85_target_structure.tsv'
	feature_df = pd.read_csv(feature_file, sep = '\t', header = None)	
	# iterate by row of predictive features 
	pred_prob_dict = {}
	for i in range(0, feature_df.shape[0]):
		# obtain the target and predictive fingerprints of current target     
		i_target = feature_df.iloc[i,0]
		i_features = feature_df.iloc[i,1].split(',')
		# load trained classifier of current target, implement the classifier to predict target binding probability 
		it_model_file = 'data/binding/' + i_target + '.joblib'
		i_classifier = joblib.load(it_model_file)
		pred_prob_dict[i_target] = i_classifier.predict_proba(structure_df[i_features])[:,1]

	## 2. Collect predicted target binding probability profile of query compounds from all targets, output in data frame 
	pred_prob_df = pd.DataFrame(pred_prob_dict)
	pred_prob_df.index = structure_df.index

	return pred_prob_df
	
