# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains data-formatting functions used for DTox model training


## Module
import numpy as np
import pandas as pd
import torch
from torch.utils.data.dataset import Dataset
from sklearn.model_selection import train_test_split


## This function formats feature-response array data into tensors used for DTox model training 
class DTox_dataformat(Dataset):
	## 0. Input arguments 
		# feature_data: array that contains input feature data 
		# label_data: array that contains input label/response data 
	
	## 1. Convert feature and label arrays into tensors 	
	def __init__(self, feature_data, label_data):
		super(DTox_dataformat, self).__init__()
		self.features = torch.tensor(feature_data, dtype = torch.float)
		labels = torch.tensor(label_data, dtype = torch.float)
		self.labels = labels.view(labels.shape[0], 1)

	## 2. Get feature and label data by index
	def __getitem__(self, index):
		feature = self.features[index]
		label = self.labels[index]
		return feature, label
	
	## 3. Obtain number of data samples  
	def __len__(self):
		return len(self.labels)


## This function converts training, testing data into the format used for DTox model training 
def format_dtox_combine_data(combine_df, label_col, N_batch = 32, test_prop = 0.125):
	## 0. Input arguments
		# combine_df: data frame tat contains combined training and testing data 
		# label_col: name of column that contains label/response data
		# N_batch: number of mini-batches 
		# test_prop: proportion of testing samples among combined training and testing data

	## 1. Format training and testing data  
	# split combined training and testing data into training and testing according to test_prop
	X_combine, y_combine = combine_df.drop(label_col, axis = 1).values, combine_df[label_col].values
	X_train, X_test, y_train, y_test = train_test_split(X_combine, y_combine, test_size = test_prop, random_state = 0, stratify = y_combine)
	# format feature and label data of training data, then generate data loader accroding to N_batch
	dtox_train_data = DTox_dataformat(X_train, y_train)
	dtox_train_data_loader = torch.utils.data.DataLoader(dtox_train_data, batch_size = N_batch, shuffle = True)
	# format feature and label data of testing data  
	dtox_test_data = DTox_dataformat(X_test, y_test)
	
	return dtox_train_data_loader, dtox_test_data


## This function converts validation data into the format used for DTox model training  
def format_dtox_valid_data(valid_df, label_col):
	## 0. Input arguments
		# valid_df: data frame that contains validation data
		# label_col: name of column that contains label/response data

	## 1. Format validation data 
	# format feature and label data of validation data 
	X_valid, y_valid = valid_df.drop(label_col, axis = 1).values, valid_df[label_col].values
	dtox_valid_data = DTox_dataformat(X_valid, y_valid) 

	return dtox_valid_data


