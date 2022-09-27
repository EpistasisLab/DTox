# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script learns and evaluates DTox model


## Module
import sys
import torch
import numpy as np
import pandas as pd
sys.path.insert(0, 'code/')
import dtox_data
import dtox_learning


## This function implements DTox model training  
def dtox(dtox_combine_df, label_col_name, out_folder, root_process = 'GE+IS+M+ST', min_pathway_size = 5, max_module_size = 20, auxiliary_alpha = 0.5, l2_lambda = 0.0001):
	## 0. Input arguments 
		# dtox_combine_df: data frame that contains combined training and testing data  
		# label_col_name: name of column that contains label/response data
		# out_folder: name of folder to store output files 
		# root_process: abbreviated name of root biological processes (see row 'root_name_annotations' of 'data/reactome/root_file_map.tsv' for values)  
		# min_pathway_size: minimal size of pathways included in DTox model 
		# max_module_size: maximal size of node modules in DTox neural network model
		# auxiliary_alpha: coefficient for auxiliary loss in loss function (see 'dtox_loss')
		# l2_lambda: coefficient for L2 regularization 

	## 1. Prepare for model training  
	# format input trainig, testing, and validation data 
	training_data_loader, testing_data = dtox_data.format_dtox_combine_data(dtox_combine_df, label_col_name)
	# specify name for output model file 
	output_model_name = out_folder + 'model.pt'
	# check whether GPU training is avaiable, if not use CPU training
	device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
	# specify name for the input ontology hierarchy files 
	root_map_file = 'data/reactome/root_file_map.tsv'
	root_map_df = pd.read_csv(root_map_file, sep = '\t', header = 0)	 
	root_id = root_map_df[root_map_df.root_name_annotations == root_process].file_id.values[0]
	rt_root_file = 'data/reactome/hierarchy/' + root_id + '_ps_5_re_0_st_0_root.tsv'
	rt_connection_file = 'data/reactome/hierarchy/' + root_id + '_ps_5_re_0_st_0_knowledge_by_node.tsv'
	rt_size_file = 'data/reactome/hierarchy/' + root_id + '_ps_5_re_0_st_0_node_size.tsv'
	rt_layer_file = 'data/reactome/hierarchy/' + root_id + '_ps_5_re_0_st_0_layer.tsv'

	## 2. learn DTox model with training and testing data (training data for learning parameters, testing data for implementing early stop), save model to output model file 
	torch.manual_seed(0)
	hierarchy_info, trained_model, loss, training_summary_df = dtox_learning.train_dtox_model(rt_root_file, rt_connection_file, rt_size_file, rt_layer_file, min_pathway_size, max_module_size, training_data_loader, testing_data, auxiliary_alpha, l2_lambda, dtox_model_name = output_model_name, dtox_device = device)

	return hierarchy_info, trained_model, loss, training_summary_df 


## This function implements DTox model evaluation 
def dtox_eval(dtox_combine_df, dtox_valid_df, label_col_name, trained_model, loss): 	
	## 0. Input arguments 
		# dtox_combine_df: data frame that contains combined training and testing data  
		# dtox_valid_df: data frame that contains validation data 
		# label_col_name: name of column that contains label/response data
		# trained_model: trained DTox model
		# loss: DTox loss function 

	## 1. Prepare for model evaluation 
	# format input learning and validation data  
	combined_data = dtox_data.format_dtox_valid_data(dtox_combine_df, label_col_name)
	validation_data = dtox_data.format_dtox_valid_data(dtox_valid_df, label_col_name)
	# check whether GPU training is avaiable, if not use CPU training 	
	device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

	## 2. Evaluate DTox model
	# evaluate learned DTox model on combined training and testing data (training performance) 
	combined_perf = dtox_learning.evaluate_dtox_model(trained_model, loss, combined_data, dtox_device = device)
	# evaluate learned DTox model on validation data (validation performance)
	validation_perf = dtox_learning.evaluate_dtox_model(trained_model, loss, validation_data, dtox_device = device)

	return combined_perf, validation_perf


