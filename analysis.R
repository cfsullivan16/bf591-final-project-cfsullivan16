#!/usr/bin/Rscript
## Caroline Sullivan
## BF591 Final Project
## Analysis

libs <- c("GEOquery", "tidyverse", "ggVennDiagram", "BiocManager",
          "DESeq2")

for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}

get_data <- function() {
  # get the metadata
  gse3 <- GEOquery::getGEO("GSE150450", GSEMatrix = TRUE)
  meta3 <- pData(gse3[[1]])
  
  # divide samples up to prep for DESEQ2 analysis
  male_titles <- meta3[meta3$'Sex:ch1' == 'male','title'] 
  # remove the male outlier
  male_titles <- male_titles[! male_titles %in% c('CAL2M')]
  female_titles <- meta3[meta3$'Sex:ch1' == 'female','title']
  larvae_titles <- meta3[meta3$'lifestage:ch1' == 'larvae','title']
  
  # load in counts
  counts <- read.csv('data/GSE150450_gene_count_matrix.csv')
  # turn gene_ids into rownames
  rownames(counts) <- counts$gene_id
  counts = subset(counts, select = -c(gene_id))
  
  male_counts <- counts[,male_titles]
  female_counts <- counts[,female_titles]
  larvae_counts <- counts[,larvae_titles]
  
  return(list(male_counts, female_counts, larvae_counts))
}

get_md <- function() {
  meta <- read.table('data/meta.csv', sep=',', header=TRUE, row.names = 1)
  return(meta)
}

##### DESEQ2 #####

# get the col_data to input into the DESeq2 analysis
get_col_data <- function(counts, md) {
  # input column name data from the counts df and experiment type
  col_data <- data.frame(samplename = colnames(counts))
  
  # add rownames
  row.names(col_data) <- colnames(counts)
  
  # add condition info
  col_data[c('condition')] <- (md %>%
    dplyr::filter(title %in% colnames(counts)))$'treatment.ch1'
  
  # add timepoint info
  col_data[c('timepoint')] <- (md %>%
    dplyr::filter(title %in% colnames(counts)))$'timepoint.ch1'
  
  # make sure columns are factors
  col_data$timepoint <- as.factor(col_data$timepoint)
  col_data$condition <- as.factor(col_data$condition)
  
  return(col_data[c('condition', 'timepoint')])
}

# run DESeq2
run_deseq <- function(count_dataframe, coldata, count_filter, condition_name) {
  # make deseqdataset object
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_dataframe),
                                colData = coldata,
                                design = ~ condition)
  
  # filter for rows with at least the number of reads indicated by count_filter
  keep <- rowSums(counts(dds)) >= count_filter
  dds <- dds[keep,]
  
  # account for timepoint
  design(dds) <- formula(~ timepoint + condition + timepoint:condition)
  
  # do differential expression analysis
  dds <- DESeq(dds, test="LRT", reduced=~ timepoint + condition)
  
  # save results
  results <- results(dds, name=condition_name)

  return(results)
}

# get genes that are significantly differentially expressed (padj < 0.01)
de_genes <- function(deseq_results){
  deseq_results <- deseq_results[!is.na(deseq_results$padj),]
  de_genes <- row.names(deseq_results[deseq_results$padj < 0.01,])
  return(de_genes)
}

