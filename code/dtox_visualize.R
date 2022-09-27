# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script uses visNetwork package to visualize the flow of relevance along VNN paths of interest between query compound, hidden pathway modules, and the outcome of interest  


## functions
library(visNetwork);
library(RColorBrewer);


## 0. Input arguments
Args			<- commandArgs(T);
query_compound		<- Args[1];		# name of query compound 
display_mode		<- Args[2];		# name of display mode ('simple': only display the paths identified by DTox, 'complex': display all paths that involve proteins connected to the query compound)
out_folder		<- Args[3];		# name of folder to store output files 
uniprot_file		<- "data/reactome/uniprot_id_map.tsv";
reactome_file		<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/reactome/UniProt2Reactome_All_Levels.txt";

## 1. Obtain names of input files  
# obtain name of input DTox module index file  
all_file_names <- readLines(paste(out_folder, "file_names.txt", sep = ""));
module_id_file <- strsplit(all_file_names[[1]], "node_file:")[[1]][[2]];
# obtain name of input query DTox paths file  
all_path_file <- strsplit(all_file_names[[2]], "path_file:")[[1]][[2]];
# obtain name of input query compound DTox module relevance score file 
compound_module_file <- strsplit(all_file_names[[3]], "pathway_relevance_file:")[[1]][[2]];
# obtain name of input query compound significant DTox path file 
compound_path_file <-  strsplit(all_file_names[[4]], "relevance_fdr_file:")[[1]][[2]];

## 2. Process input query compound info and DTox model interpretation result  
# read in significant DTox paths of compounds for DTox model as data frame, select rows of interest to obtain significant DTox paths identified for the query compound  
compound_path_df <- read.delim(file = compound_path_file, header = T, sep = "\t");
cpd_id <- which(compound_path_df$cid %in% query_compound);
query_path <- as.character(compound_path_df$path_id[cpd_id]); 
# obtain IDs of proteins linked to the query compound  
query_proteins <- sapply(query_path, function(qp){
	qp_s <- strsplit(qp, "_")[[1]];
	qs_len <- length(qp_s);
	return(qp_s[[qs_len]]);
});
query_proteins <- unique(query_proteins);
# 
uniprot_df <- read.delim(file = uniprot_file, header = T, sep = "\t");
qp_id <- which(uniprot_df$protein_uniport %in% query_proteins);
query_protein_udf <- uniprot_df[qp_id, ];
colnames(query_protein_udf) <- c("node_id", "node_name");
# read in all DTox paths, iterate by path to check overlap with query proteins and obtain paths of interest for plotting
if(display_mode == "simple")	all_path_overlap <- query_path
if(display_mode == 'complex'){
	all_paths <- readLines(all_path_file)
	all_path_overlap <- lapply(all_paths, function(ap){
		# separate path name by underscore to obtain the pathways along the current path, check whether query proteins are among them 
		ap_s <- strsplit(ap, "_")[[1]]
		ap_inter <- intersect(ap_s, query_proteins)
		if(length(ap_inter) == 0)	return(logical(0))
		else return(ap)
	})
}
all_path_overlap <- unlist(all_path_overlap);
# obtain the IDs of all pathways involved in paths of interest  
apo_elements <- lapply(all_path_overlap, function(apo) strsplit(apo, "_")[[1]]);
apo_elements <- unique(unlist(apo_elements)); 
# read in pathway annotation from Reactome 
reactome_df <- read.delim(file = reactome_file, header = F, sep = "\t");
reactome_df <- unique(reactome_df[ ,c("V2", "V4")]);
# obtain the IDs of pathways of interest  
ae_id <- which(reactome_df$V2 %in% apo_elements);
apo_df <- reactome_df[ae_id, ];
colnames(apo_df) <- c("node_id", "node_name");
network_node_df <- rbind(data.frame("node_id" = query_compound, "node_name" = query_compound), query_protein_udf, apo_df, data.frame("node_id" = "root", "node_name" = "Outcome"));
# read in index of DTox module as data frame, obtain the DTox model index of query proteins and pathways    
module_id_df <- read.delim(file = module_id_file, header = T, sep = "\t");
network_node_module <- sapply(network_node_df$node_id, function(nndni){
	nndni_id <- which(module_id_df$node_name %in% nndni);
	if(length(nndni_id) == 0)	return(NA)
	else	return(module_id_df$node[[nndni_id]])
});
# read in DTox module relevance scores of compounds for DTox model as data frame, select row/columns of interest to obtain module relevance scores of query compound with respect to query proteins and pathways  
compound_module_df <- read.delim(file = compound_module_file, header = T, sep = "\t");
qc_id <- which(compound_module_df$X %in% query_compound)[[1]];
network_node_score <- sapply(network_node_module, function(nnm){
	if(is.na(nnm) == T)     nnm_score <- 0
	else{
		nnm_score <- compound_module_df[qc_id, nnm+2]
		# convert negative scores to 0 (as in our interpretation method)
		if(nnm_score < 0)	nnm_score <- 0
	}
	return(nnm_score);
});

## 2. Tune network edge parameters for plotting  
# iterate by DTox paths of interest, break each path into network edges for plotting   
all_path_edges <- lapply(all_path_overlap, function(apo){
	# add root and compound nodes to the current path 
	apo_s <- strsplit(apo, "_")[[1]];
	apo_s <- c("root", apo_s, query_compound);
	# break the path into edges that connect adjacent entities 
	as_len <- length(apo_s);
	apo_edges <- mapply(function(al) apo_s[al:(al-1)], as_len:2);
	return(t(apo_edges));
});
# iterate by DTox paths of interest, check whether each path was identified with the query compound     
all_path_query <- sapply(all_path_overlap, function(apo){
	qp_apo <- sapply(query_path, function(qp) length(grep(apo, qp)));
	return(sum(qp_apo));
});
# decide whether to draw each DTox path of interst in dashed or solid lines
all_edge_dash <- mapply(function(apq, ape){
	# if path was not identified with the query compound, draw in dashed lines, otherwise draw in solid lines 
	apq_l <- (apq == 0)
	apeq <- rep(apq_l, nrow(ape));
	return(apeq);
}, all_path_query, all_path_edges);
# aggregate edges from all DTox paths of interest, obtain unique network edge set 
all_path_edges <- do.call(rbind, all_path_edges)
edges1 <- data.frame(from = all_path_edges[,1], to = all_path_edges[,2], dashes = unlist(all_edge_dash));
edges1$from_to <- mapply(function(ef, et) paste(ef, et, sep = ","), edges1$from, edges1$to);
unique_eft <- unique(edges1$from_to);
# specify whether to draw each edge in dashed or solid line 
ue_edges <- mapply(function(ue){
	# if all DTox paths associated with the edge were not identified with the query compound, draw in dashed lines, otherwise draw in solid lines 
	ue_id <- which(edges1$from_to %in% ue);
	ue_dash <- prod(edges1$dashes[ue_id]);
	return(c(ue_id[[1]], ue_dash));
}, unique_eft);
# specify edge color  
edge_col <- sapply(ue_edges[2, ], function(ue2){
	# if all DTox paths associated with the edge were not identified with the query compound, draw in grey, otherwise draw in purple 
	if(ue2 > 0)	return('grey')
	else	return('purple')
});
# store edge connection, shape, and color info in data frame 
edges <- data.frame(from = all_path_edges[ue_edges[1,], 1], to = all_path_edges[ue_edges[1,], 2], dashes = ue_edges[2, ] > 0, color = edge_col);

## 3. Tune network node parameters for plotting 
# adjust node label name to accommodate space between network layers  
network_node_name <- sapply(as.character(network_node_df$node_name), function(nndnn){
	nndnn_s <- strsplit(nndnn, " ")[[1]];
	ns_len <- length(nndnn_s);
	if(ns_len == 1)	return(nndnn)
	else if(ns_len < 6){
		cut_len <- ceiling(ns_len/2)
		nndnn1 <- paste(paste(nndnn_s[1:cut_len], collapse = " "), paste(nndnn_s[(cut_len+1):ns_len], collapse = " "), sep = "\n")
		return(nndnn1)
	}
	else{
		cut_len1 <- ceiling(ns_len/3)
		cut_len2 <- cut_len1*2
		nndnn2 <- paste(paste(nndnn_s[1:cut_len1], collapse = " "), paste(nndnn_s[(cut_len1+1):cut_len2], collapse = " "), paste(nndnn_s[(cut_len2+1):ns_len], collapse = " "), sep = "\n")		
	}
});
# define node color based on the scale of pathway relevance scores
nns_min <- min(network_node_score[which(network_node_score > 0)]);
network_node_score[which(network_node_score == 0)] <- nns_min/10;
network_node_score <- log10(network_node_score);
node_pal <- colorRampPalette(c('white', 'purple', 'purple1'));
network_node_col <- node_pal(10)[as.numeric(cut(network_node_score, breaks = 10))];
# define node shape based on node class: triangle for query compound, circle for pathways, square for root outcome  
node_shape <- sapply(as.character(network_node_df$node_name), function(nndnn){
	if(nndnn == query_compound)	return("triangle")
	else if(nndnn == "Outcome")	return("box")
	else if(nndnn %in% as.character(query_protein_udf$node_name))	return("diamond")
	else return("dot")
});
# define node border color, grey for pathways, black for other nodes  
network_node_border <- rep("grey", length(network_node_df$node_id));
names(network_node_border) <- network_node_df$node_id;
network_node_border[[query_compound]] <- "black";
network_node_border[["root"]] <- "black"
# store node label, color, shape, and border color info in data frame  
nodes <- data.frame(id = network_node_df$node_id, label = network_node_name, color.background = network_node_col, color.border = network_node_border, shape = node_shape);

## 3. Make the network plot based on specified parameters, showing the flow of relevance along nodes/edges in the network 
network_plot_file <- paste(out_folder, query_compound, ".html", sep = "");
network <- visNetwork(nodes, edges) %>%
	visNodes(size = 15, font = '8px helvetica black') %>%
	visEdges(arrows = "to") %>%
	visHierarchicalLayout(direction = "LR", sortMethod = 'directed', levelSeparation = 180) %>%
	visEvents(stabilizationIterationsDone = "function(){this.setOptions({physics:false});}") %>%
	visExport(type = "pdf", name = network_plot_file)
visSave(network, file = network_plot_file);
