## Function for preparing gene lists for gsea or over representation analysis
## Developed by George William Kagugube 
## Date: 29 Jan 2024 
## Usage, read teh file into R preferably using source (file_name)
## prepare an object from DESeq of the differentially expressed genes
## This is given as an  data frame to the function here

# Load the libraries that will be needed for the analysis herein 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(clusterProfiler) # for PEA analysis
library(org.Dr.eg.db)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(GOSemSim)
library(DOSE)
library(pathview)
library(enrichplot)
library(AnnotationDbi)
library(AnnotationHub)

## The functions start from here
transcriptome_annotation <- function(data_obj){
  # lOAD THE REQUIRED LIBRARIES HERE
  library("AnnotationDbi")
  library("org.Dr.eg.db")
  
  ### Add to the data structure of the gene symbols and entrezid
  data_obj$ENTREZID <- mapIds(org.Dr.eg.db,
                              keys = row.names(brain_mutant_evu),
                              column = c("ENTREZID"),
                              keytype = "ENSEMBL",
                              multiVals = "first")
  
  data_obj$Symbols <- mapIds(org.Dr.eg.db,
                             keys = row.names(brain_mutant_evu),
                             column = c("SYMBOL"),
                             keytype = "ENSEMBL",
                             multiVals = "first")
  
  ## Return the annotated data object here 
  return(data_obj)
}

# Function for generating the 
enrich_ont <- function(dge_file, method="gsea"){
  
  if (method == "gsea"){
    # Display a message of what analysis is being conducted here
    print ("You're analysing your gene list using gsea")
    # Select the fold change from the 
    original_gene_list <- dge_file$log2FoldChange
    names(original_gene_list) <- rownames(dge_file)
    original_gene_list <- na.omit(original_gene_list)
    gene_list <- sort(original_gene_list, decreasing = TRUE)
    
    return(gene_list)
  }else{
    print ("You're analysing your gene list using over representation")
    # we want the log2 fold change 
    
    gene_list <- dge_file[(dge_file$baseMean > 50 & dge_file$padj < 0.05),]
    genes <- rownames(gene_list)
    
    return(genes)
  }
}

### Pathway analysis gene list preperation 
kegg_pathway_analysis <- function(dge_file, method = "gsea"){
  
  # Select the fold change from the 
  original_gene_list <- dge_file$log2FoldChange
  names(original_gene_list) <- rownames(dge_file)
  
  if (method == "gsea"){
    print ("You're analysing your gene list using gsea")
    ## Convert the ensembl IDs to kegg API naming convention 
    ids <- bitr(geneID = names(original_gene_list),
                     fromType = "ENSEMBL",
                     toType = "ENTREZID",
                     OrgDb = "org.Dr.eg.db",
                     drop = TRUE)
    
    # Remove any duplicates from the list above 
    dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]),]
    
    # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
    df2 = dge_file[rownames(dge_file) %in% dedup_ids$ENSEMBL,]
    
    # Create a new column in df2 with the corresponding ENTREZ IDs
    rownames(dedup_ids) <- dedup_ids$ENSEMBL 
    
    df2 <- merge(dedup_ids, dge_file, by="row.names")
    
    ## Select the desired genes here 
    kegg_gene_list <- df2$log2FoldChange
    names(kegg_gene_list) <- df2$ENTREZID
    kegg_gene_list <- na.omit(kegg_gene_list)
    kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
    
    kegg_gene_list <- kegg_gene_list[!duplicated(kegg_gene_list)]
    
    return(kegg_gene_list)
    
  }else{
    print ("You're analysing your gene list using over representation")
    ids<-bitr(names(original_gene_list), 
              fromType = "ENSEMBL", 
              toType = "ENTREZID", 
              OrgDb="org.Dr.eg.db")
    
    dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
    df2 = dge_file[rownames(dge_file) %in% dedup_ids$ENSEMBL,]
    # Create a new column in df2 with the corresponding ENTREZ IDs
    df2$ENTREZID = dedup_ids$ENTREZID
    
    # Create a vector of the gene unuiverse
    kegg_gene_list <- df2$log2FoldChange
    
    # Name vector with ENTREZ ids
    names(kegg_gene_list) <- df2$ENTREZID
    
    # omit any NA values 
    kegg_gene_list<-na.omit(kegg_gene_list)
    
    # sort the list in decreasing order (required for clusterProfiler)
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
    
    # Exctract significant results from df2
    kegg_sig_genes_df = subset(df2, padj < 0.05)
    
    # From significant results, we want to filter on log2fold change
    kegg_genes <- kegg_sig_genes_df$log2FoldChange
    
    # Name the vector with the CONVERTED ID!
    names(kegg_genes) <- kegg_sig_genes_df$ENTREZID
    
    # omit NA values
    kegg_genes <- na.omit(kegg_genes)
    
    # filter on log2fold change (PARAMETER)
    kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 0.5]
    
    return(kegg_genes)
  }
}

##### Functional analysis by calling gse
gseFunc <- function(gene_list, item="ALL"){
  gse <- gseGO(geneList = gse_gene_list,
               ont = item,
               keyType = "ENSEMBL",
               minGSSize = 5,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               eps = 0,
               OrgDb = "org.Dr.eg.db",
               pAdjustMethod = "fdr")
  return(gse)
}

oraFunc <- function (gene_list, item = "ALL") {
  go_enrich <- enrichGO(gene = ora_gene_list,
                        OrgDb = "org.Dr.eg.db", 
                        keyType = 'ENSEMBL',
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = 0.10)
  return(go_enrich)
}

### Visualisation functions
plotting_pathways <- function(obj, mode = "gsea") {
  
  if (mode == "gsea"){
    plot <- dotplot(obj, 
                    color = "p.adjust",
                    showCategory=10, 
                    font.size = 20,
                    split=".sign") + facet_grid(.~.sign)
    
    return(plot)
  }else{
    plot <- dotplot(obj, 
                    color = "p.adjust",
                    showCategory=10, 
                    font.size = 20)
    
    return(plot)
  }
  plot <- dotplot(obj, 
                  color = "p.adjust",
                  showCategory=10, 
                  font.size = 17,
                  split=".sign") + facet_grid(.~.sign)
  
  return(plot)
}

emap_plotting <- function(obj){
  ego2 <- pairwise_termsim(obj)
  plot <- emapplot(ego2)
  
  return(plot)
}