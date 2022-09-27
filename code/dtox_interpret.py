# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script contains implements layer-wise relevance propagation to evaluate relevance of DTox paths


## Module
import os
import sys
import numpy as np
import pandas as pd
import torch
sys.path.insert(0, 'code/')
import dtox_hierarchy
import dtox_nn
import dtox_data
import dtox_lrp
import dtox


def dtox_interpret(query_data_df, dtox_model, dtox_combine_df, label_col_name, out_folder, prop_rule = 'gamma-epsilon', rule_factor1 = 0.001, rule_factor2 = 0.1, N_null_models = 200, root_process = 'GE+IS+M+ST', min_pathway_size = 5, max_module_size = 20, auxiliary_alpha = 0.5, l2_lambda = 0.0001):
	## 0. Input arguments 
		# query_data_df: data frame that contains the data instances for interpretation 
		# dtox_model: trained DTox model
		# dtox_combine_df: data frame that contains combined training and testing data  
		# label_col_name: name of column that contains label/response data
		# out_folder: name of folder to store output files  
		# prop_rule: propagation rule to be implemented ('gamma-epsilon' for generic rule, 'alpha-beta' for αβ rule)   
		# rule_factor1: first factor in LRP rules ('gamma' in generic rule, 'alpha' in αβ rule)  
		# rule_factor2: second factor in LRP rules ('epsilon' in generic rule, 'beta' in αβ rule)  
		# N_null_models: number of sampling times to compute the empirical P-value of path relevance scores   
		# root_process: abbreviated name of root biological processes (see row 'root_name_annotations' of 'data/reactome/root_file_map.tsv' for values)  
		# min_pathway_size: minimal size of pathways included in DTox model 
		# max_module_size: maximal size of node modules in DTox neural network model
		# auxiliary_alpha: coefficient for auxiliary loss in loss function (see 'dtox_loss')
		# l2_lambda: coefficient for L2 regularization 
	
	## 1. Prepare for model interpretation  
	# format featuer and label data of query dataset 
	X_query, y_query = query_data_df.drop(label_col_name, axis = 1).values, query_data_df[label_col_name].values
	query_data = dtox_data.DTox_dataformat(X_query, y_query)
	query_feature_data = query_data[: , :][0]
	# obtain name of node information file 
	root_map_file = 'data/reactome/root_file_map.tsv'
	root_map_df = pd.read_csv(root_map_file, sep = '\t', header = 0)
	root_id = root_map_df[root_map_df.root_name_annotations == root_process].file_id.values[0]
	rt_node_file = 'data/reactome/hierarchy/' + root_id + '_ps_5_re_0_st_0_node.tsv'
	
	## 2. Implement layer-wise relevance propagation (LRP) process to compute the relevance of all neurons/modules in learned DTox model 
	hidden_relevance_df, pathway_relevance_df = dtox_lrp.lrp(dtox_model, query_feature_data, prop_rule, float(rule_factor1), float(rule_factor2), 0, 1)
	hidden_relevance_df.index = query_data_df.index
	pathway_relevance_df.index = query_data_df.index
	
	## 3. Compute observed relevance scores for all DTox paths connecting input layer to root module  
	path_relevance_df = dtox_lrp.compute_path_relevance_score(dtox_model, pathway_relevance_df)
	# read in data frame that contains node name info 
	node_map_df = pd.read_csv(rt_node_file, sep = '\t', header = 0)
	# iterate by path, name each path by names of nodes in the path  	
	path_name_dict = {}
	for prdc in path_relevance_df.columns:
		# link node numbers/names of current path by underscores respectively, store in dictionary 
		prdc_node_id = [int(pr) for pr in prdc.split('_')]
		prdc_path_name = '_'.join(node_map_df.node_name[prdc_node_id].values)
		path_name_dict[prdc] = prdc_path_name
	
	## 4. Compute empirical p-values of module relevance scores for all DTox path-compound pairs 
	# iterate by trained null DTox model
	null_relevance_list = []
	null_file_folder = out_folder + 'null/'
	os.mkdir(null_file_folder)
	for dnf in range(0, N_null_models):
		# shuffle the outcome labels of samples among training data
		dnf_combine_df = dtox_combine_df.copy()
		dnf_combine_df[label_col_name] = dnf_combine_df[label_col_name].sample(frac = 1, random_state = dnf).values
		# define structure of whole neural network of null model 
		dnf_model_file = null_file_folder + 'null_' + str(dnf + 1) + '_'
		_, dnf_model, _, _ = dtox.dtox(dnf_combine_df, label_col_name, out_folder = dnf_model_file, root_process = root_process, min_pathway_size = min_pathway_size, max_module_size = max_module_size, auxiliary_alpha = auxiliary_alpha, l2_lambda = l2_lambda)
		# implement LRP to compute the relevance of all modules in null DTox model 
		_, dnf_pathway_relevance_df = dtox_lrp.lrp(dnf_model, query_feature_data, prop_rule, float(rule_factor1), float(rule_factor2), 0, 1)
		dnf_pathway_relevance_df.index = query_data_df.index
		# compute null relevance scores for all DTox paths 
		dnf_path_relevance_df = dtox_lrp.compute_path_relevance_score(dnf_model, dnf_pathway_relevance_df)
		null_relevance_list.append(dnf_path_relevance_df)
	# compute empirical p-value of relevance scores for all DTox path-compound pairs
	relevance_fdr_df = dtox_lrp.compute_path_relevance_pvalue(path_relevance_df, null_relevance_list) 
	# iterate by path, substitute node number with node name   
	relevance_path_names = []
	for rpv in relevance_fdr_df.path_id.values:
		relevance_path_names.append(path_name_dict[rpv])
	relevance_fdr_df.path_id = relevance_path_names

	## 5. Output result files 
	# write all named paths to output txt file
	path_file = out_folder + 'all_paths.txt'
	path_file_o = open(path_file, 'w')
	for lpv in list(path_name_dict.values()):
		path_file_o.write('%s\n' % lpv)
	path_file_o.close()
	# write data frame of node module relevance scores of compounds to output file 
	pathway_relevance_file = out_folder + 'module_relevance.tsv'
	pathway_relevance_df.columns = node_map_df.node_name.values 
	pathway_relevance_df.to_csv(pathway_relevance_file, sep = '\t')
	# write data frame that contains p-values of DTox path-compound pairs to output tsv file 
	relevance_fdr_file = out_folder + 'path_relevance_pv.tsv'
	relevance_fdr_df.to_csv(relevance_fdr_file, sep = '\t', float_format = '%.5f', index = False)
	# write file names to file
	file_name = ['node_file:' + rt_node_file, 'path_file:' + path_file, 'pathway_relevance_file:' + pathway_relevance_file, 'relevance_fdr_file:' + relevance_fdr_file] 
	name_file_o = open(out_folder + 'file_names.txt', 'w')
	for fn in file_name:
		name_file_o.write('%s\n' % fn)
	name_file_o.close()

	return pathway_relevance_df, relevance_fdr_df


