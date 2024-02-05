gene_list <- function(dge_file, method){
	# Extract the gene list for gsea analysis
	if (method == "gsea"){
	# Display a message of what analysis is being conducted here
	print ("You're analysing your gene list using gsea")

	original_gene_list <- dge_file$log2FoldChange
	names(original_gene_list) <- rownames(original_gene_list)

	# omit any NA values 
	gene_list<-na.omit(original_gene_list)

	# sort the list in decreasing order (required for clusterProfiler)
	gene_list = sort(gene_list, decreasing = TRUE)

	return(gene_list)
	}else{
		print ("You're analysing your gene list using over representation")
	# we want the log2 fold change 
	original_gene_list <- dge_file$log2FoldChange

	# name the vector
	names(original_gene_list) <- rownames(dge_file)

	# omit any NA values 
	gene_list<-na.omit(original_gene_list)

	# sort the list in decreasing order (required for clusterProfiler)
	gene_list = sort(gene_list, decreasing = TRUE)

	# Exctract significant results (padj < 0.05)
	sig_genes_df = subset(dge_file, padj < 0.05)

	# From significant results, we want to filter on log2fold change
	genes <- sig_genes_df$log2FoldChange

	# Name the vector
	names(genes) <- rownames(sig_genes_df)

	# omit NA values
	genes <- na.omit(genes)

	# filter on min log2fold change (log2FoldChange > 2), 
	# this can be changed accordingly, for my data I choose a cut off of 0.2
	genes <- names(genes)[abs(genes) > 0.2]

	return (genes)

	}
}

pathway_analysis <- function(dge_file, method){
	# Check for the method chosen by the scientist here 
	if (method == "gsea"){

		# Extract the genes of interest here. This is the entire list of genes from 
		original_gene_list <- dge_file$log2FoldChange
	    names(original_gene_list) <- rownames(original_gene_list)


		# Convert gene IDs for gseKEGG function
		# We will lose some genes here because not all IDs will be converted
		ids<-bitr(names(original_gene_list), 
			fromType = "ENSEMBL", 
			toType = "ENTREZID", 
			OrgDb=organism)
		 # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
		dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

		# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
		df2 = df[df$X %in% dedup_ids$ENSEMBL,]

		# Create a new column in df2 with the corresponding ENTREZ IDs
		df2$Y = dedup_ids$ENTREZID

		# Create a vector of the gene unuiverse
		kegg_gene_list <- df2$log2FoldChange

		# Name vector with ENTREZ ids
		names(kegg_gene_list) <- df2$Y

		# omit any NA values 
		kegg_gene_list<-na.omit(kegg_gene_list)

		# sort the list in decreasing order (required for clusterProfiler)
		kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

		return (kegg_gene_list)

		}else{

			# Convert gene IDs for enrichKEGG function
			# We will lose some genes here because not all IDs will be converted
			ids<-bitr(names(original_gene_list), 
				fromType = "ENSEMBL", 
				toType = "ENTREZID", 
				OrgDb="org.Dr.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
			dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

			# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
			df2 = df[df$X %in% dedup_ids$ENSEMBL,]

			# Create a new column in df2 with the corresponding ENTREZ IDs
			df2$Y = dedup_ids$ENTREZID

			# Create a vector of the gene unuiverse
			kegg_gene_list <- df2$log2FoldChange

			# Name vector with ENTREZ ids
			names(kegg_gene_list) <- df2$Y

			# omit any NA values 
			kegg_gene_list<-na.omit(kegg_gene_list)

			# sort the list in decreasing order (required for clusterProfiler)
			kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

			# Exctract significant results from df2
			kegg_sig_genes_df = subset(df2, padj < 0.05)

			# From significant results, we want to filter on log2fold change
			kegg_genes <- kegg_sig_genes_df$log2FoldChange

			# Name the vector with the CONVERTED ID!
			names(kegg_genes) <- kegg_sig_genes_df$Y

			# omit NA values
			kegg_genes <- na.omit(kegg_genes)

			# filter on log2fold change (PARAMETER)
			kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]

			return (kegg_genes)

		}
}