################################################################################
##########  Set up the environment for the analysis ############################
################################################################################
# Clear the environmental variables for a new analysis here
rm(list = ls())

# Set the working directory here
setwd("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data")

# Check that your are in the right environment
getwd()

# Check the contents of the working environment here
dir()
################################################################################
### Load the required files here
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
################################################################################
####### lOAD THE DAT COUNT MATRIX AND META/SAMPLE INFORMATION ##################
################################################################################
################# Load the brain data here #####################################
countData <- read.csv("countData_brain.csv",
                      sep = ',', header = TRUE,row.names = 1)
sample_info <- read.csv("brain_sample_information.csv",
                        sep = ',', header = TRUE,row.names = 1)
annotation <- read.csv("annotation.csv",
                       row.names = 1)

############# load Eye data information ########################################
# countData <- read.csv("countData_eyes.csv", sep = ',', header = TRUE,
#                     row.names = 1)
# sample_info <- read.csv("eye_sample_information.csv", sep = ',', header = TRUE,
#                      row.names = 1)

########### Explore the countdata and sample information here ##################
head(countData)
tail(countData)
dim(countData)
str(countData)

# Check the data structure in 
head(sample_info, 13)
tail(sample_info)
dim(sample_info)
str(sample_info)

head(annotation, 10)
######## Rearrange the columns of the countData matrix to match the rownames of 
######## the sample information table ##########################################
## Dropping some of the samples from the dataframe here
sample_info <- sample_info[c(2,3,4,5,6,8,9,10,11),]

sample_info <- sample_info[c(2,3,4,5,6,8,9,10),]
# Create a vector of the names and order of samples from the sample info
info <- c(rownames(sample_info))

# Use the information above to the selection here
countData <- countData %>%
  dplyr::select(all_of(info))

countData <- countData %>%
  dplyr::select(1:9)

###############################################################################

## Check that the columns match the rows (countData vs sample information)
names(countData)
attach(sample_info)
attach(countData)

# Check whether column names match the row name here
all(rownames(sample_info) %in% colnames(countData))
all(rownames(sample_info) == colnames(countData))

################################################################################
###### Create the DESeq object here 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sample_info,
                              design = ~ group)

dim(dds)
## Pre-filtering (not required but recommended)
## Remove all the genes with less than 10 read counts within a sample
## minimum sample replicate in the data set here 
# smallestGroupSize <- 3
# keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
# dds <- dds[keep,]
# dim(dds)
################################################################################
########### Perform some normalisation  ########################################
################################################################################
## Plot the distribution of each sample here 
# jpeg("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/figures/sample_distall.jpg",
#      width = 350, height = 350)
ggplot(data = countData) +
  geom_histogram(aes(x=Sample2), stat = "bin", bins = 200) +
  # zoom in to see the distribution more clearly
  xlim(-2, 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 10)
  )

# dev.off()

# Plot of varience vs mean to determine the appropriate distribution to model
# the data 
mean_counts <- apply(countData[, 1:length(colnames(countData))], 1, mean)
variance_counts <- apply(countData[, 1:length(colnames(countData))], 1, var)
df <- data.frame(mean_counts, variance_counts)

# jpeg("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/figures/mean_var_distall.jpg",
#      width = 350, height = 350)
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic() +
  theme(
    axis.text.y = element_text(family = "Times-Roman", size = 15),
    axis.text.x = element_text(family = "Times-Roman",size = 15),
    axis.title.y = element_text(family = "Times-Roman",size = 20),
    legend.text = element_text(family = "Times-Roman",size = 15),
    legend.title = element_text(family = "Times-Roman",size = 20)
  )
# dev.off()

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))


# Plot PCA 
# jpeg("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/figures/pcaplot1.jpg",
#      width = 350, height = 350)
z <- plotPCA(rld, intgroup="Genotype", ntop = 1000)
z + geom_label(aes(label = name)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(family = "Times-Roman", size = 15),
    axis.text.x = element_text(family = "Times-Roman",size = 15),
    axis.title.y = element_text(family = "Times-Roman",size = 20),
    legend.text = element_text(family = "Times-Roman",size = 15),
    legend.title = element_text(family = "Times-Roman",size = 20)
  )

# dev.off()

## Check with other PCAs for variation arising from the exposure vs unexposure in the groups
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(sample_info, pca$x)
head(df)

# jpeg("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/figures/pcaplotall.jpg",
#      width = 350, height = 350)
ggplot(df) + 
  geom_point(aes(x=PC1, y=PC2, color = group), size=5) +
  theme_classic() +
  theme(
    axis.text.y = element_text(family = "Times-Roman", size = 15),
    axis.text.x = element_text(family = "Times-Roman",size = 15),
    axis.title.y = element_text(family = "Times-Roman",size = 20),
    axis.title.x = element_text(family = "Times-Roman",size = 20),
    legend.text = element_text(family = "Times-Roman",size = 15),
    legend.title = element_text(family = "Times-Roman",size = 20)
  )
  
# dev.off()
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)

# jpeg("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/figures/heatmapall.jpg",
#      width = 350, height = 350)
# Plot heatmap
pheatmap(rld_cor)
# dev.off()

## Extract and the store in a csv file the normalised counts for each gene in a 
## all the samples here.
### Extract the normalised counts for use later in the functional analysis
# 
# tissue <- "brain"
# 
# # check the tissue defined and perform the right normalisation for later in the 
# # analysis here (Functional analysis)
# 
# if (tissue == "brain"){
#   dds <- estimateSizeFactors(dds)
#   sizeFactors(dds)
#   normalized_counts <- counts(dds, normalized=TRUE)
#   
#   ## Annotate the files with the normalised dataframe here
#   normalized_counts <- merge(annotation, 
#                              as.data.frame(normalized_counts),
#                              by = "row.names", all = TRUE)
#   
#   # Export the file for functional anlysis later
#   write.csv(normalized_counts, 
#             file = "/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/normalised_counts_brain_allSamples.csv")
# } else{
#   dds <- estimateSizeFactors(dds)
#   sizeFactors(dds)
#   normalized_counts <- counts(dds, normalized=TRUE)
#   
#   ## Annotate the files with the normalised dataframe here
#   normalized_counts <- merge(annotation, 
#                            as.data.frame(normalized_counts),
#                           by = "row.names", all = TRUE)
#   
#   # Export the file for functional anlysis later
#   write.csv(normalized_counts, 
#             file = "/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/normalised_counts_eyes.csv")
# }
# 
# normalized_counts
# 
# # Or we can set the comparisons as below
# dds$Treatment <- relevel(dds$Treatment, ref = "Unexposed")# Fit the statistical model here 
# dds$Genotype <- relevel(dds$Genotype, ref = "WT")
dds$group <- relevel(dds$group, ref = "WT_Unexposed")

## Run DeSeq function here 
dds <- DESeq(dds)

cbind(resultsNames(dds))[,1]

#### Check on the dispersion 
# Plot dispersion estimates
plotDispEsts(dds)

# Comparing Unexposed wild type and exposed
################################################################################
# 1. WT_unexposed vs Mut_unexposed
# 2. Mut_unexposed vs Mut_exposed
# 3. WT_unexposed vs WT_exposed
# 4. WT_unexposed vs Mut_Unexposed

################################################################################

#### Extract data with the comperisions stated above 

#res1 <- results(dds, contrast=contrast_3)
#res1 <- lfcShrink(dds, contrast = contrast_1, res = res1)
#summary(res)
res_unshrunken <- results(dds, name =  "group_Mutant_Unexposed_vs_WT_Unexposed")
summary(res_unshrunken)
res_shrucken <- lfcShrink(dds = dds, coef = 3, res = res_unshrunken)
summary(res_shrucken)
write.csv(res_shrucken, 
          file = "/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/excludedSamples/group_Mutant_Unexposed_vs_WT_Unexposed_eye_selected.csv")
plotMA(res_unshrunken)
plotMA(res_shrucken)

 ## Adjusted p values
res0rdered <- res[order(res$padj),]
res0rdered <- as.data.frame(res0rdered)
write.csv(res0rdered, 
          file = "/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/mut_unexpo_vs_mut_exposed.csv")

res2 <- results(dds, name = "group_WT_Exposed_vs_Mutant_Exposed")
summary(res2)
res2 <- lfcShrink(dds, coef = 3, res = res2)
summary(res2)
dim(res2)

res20rdered <- res2[order(res2$padj),]
res20rdered <- as.data.frame(res20rdered)
write.csv(res20rdered, 
          file = "/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/wt_expo_vs_mut_exposed.csv")


res3 <- results(dds, name =  "group_WT_Unexposed_vs_Mutant_Exposed")
summary(res3)
res3 <- lfcShrink(dds, coef = 4, res = res3)
summary(res3)

res30redered <- res3[order(res3$padj),]
res30redered <- as.data.frame(res30redered)
write.csv(res20rdered, 
          file = "/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/data/deseq_output/wt_unexpo_vs_mut_exposed.csv")

### Treatment and genotype effect 
res_treatment <- results(dds, name = "Treatment_Exposed_vs_Unexposed")
summary(res_treatment)
res_treatment <- lfcShrink(dds = dds, coef = 3, res = res_treatment)
summary(res_treatment)

res_genotype <- results(dds, name = "Genotype_Mutants_vs_WT")
summary(res_genotype)

res_interaction <- results(dds, name = "GenotypeMutants.TreatmentExposed")
summary(res_interaction)

head(coef(dds))

############ Visualisation
hist(res$pvalue)

resshrink <- lfcShrink(dds = dds, coef = 2, res = res)

shrinktable <- as.data.frame(resshrink)

shrinkTab.11 <- shrinktable

volcanoTab.11 <- shrinkTab.11 %>% 
  mutate(`-log10(pvalue)` = -log10(pvalue))

ggplot(volcanoTab.11, aes(x = log2FoldChange, y=`-log10(pvalue)`)) + 
  geom_point(aes(colour=padj < 0.05), size=1) #+
  #geom_text(data=~top_n(.x, 1, wt=-padj), aes(label=rownames))
  
vennDat <- tibble(Geneid =rownames(res)) %>%
  
  mutate(mutunexposedvsmutexposed = res$padj < 0.05 &
           !is.na(res$padj))  %>%
  mutate(wtexpovsmutexpo = res2$padj < 0.05 &
           !is.na(res2$padj)) %>%
  mutate(wtunexpovsmutexpo = res3$padj < 0.05 &
           !is.na(res3$padj))

ggvenn(vennDat, set_name_size = 3)
## Create a function of data here 
annData <- function(ann, data1) {
  ## Function takes in as input annotation and logfold datasets
  ## Remove all genes with missing padj values
  data1 <- na.omit(data1)
  
  ## Arrange the data in increasing padj values
  data1 <- data1[order(data1$padj, decreasing = FALSE),]
  
  ## Merge the annotation dataframe with the dataste from the dge expression
  merged_dge <- merge(ann, as.data.frame(data1), by = "row.names", all = TRUE)
  
  ## Return the annotated values
  return(merged_dge)
}
dim(res)
dim(annotation)

dge_list <- annData(annotation, res4)
head(dge_list, 300)
dge_list %>%
  filter(padj < 0.05) %>%
  head(150)

### Check the session information here for version controls
sessionInfo()

