#!/usr/bin/Rscript
## Caroline Sullivan
## BF591 Final Project
## Analysis

libs <- c("GEOquery", "tidyverse", "RColorBrewer", "ggVennDiagram", "BiocManager",
          "DESeq2")

for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}



################################# 
######## SET UP THE DATA ######## 
################################# 

################################
## FUNCTIONS FOR METADATA TAB ##
################################

get_md <- function() {
  # read in the metadata
  meta <- read.table('data/meta.csv', sep=',', header=TRUE, row.names = 1)
  return(meta)
}

filter_md <- function(md) {
  # filters the metadata to only the relevant columns for display in app
  keeps <- c('title', 'geo_accession', 'source_name_ch1', 'organism_ch1', 
            'molecule_ch1', 'platform_id', 'library_selection', 
            'library_source', 'library_strategy', 'lifestage.ch1', 'Sex.ch1', 
            'timepoint.ch1', 'treatment.ch1')
  
  md_filtered <- md %>%
    dplyr::select(all_of(keeps)) %>%  # get the keeps columns
    dplyr::rename('Sample name' = title,  # clean up the column names
                  'GEO Accession' = geo_accession,
                  'Source' = source_name_ch1,
                  Organism = organism_ch1,
                  Molecule = molecule_ch1,
                  'Platform ID' = platform_id,
                  'Library Selection' = library_selection,
                  'Library Source' = library_source,
                  'Library Strategy' = library_strategy,
                  'Lifestage' = 'lifestage.ch1',
                  'Sex' = 'Sex.ch1',
                  'Timepoint' = 'timepoint.ch1',
                  'Treatment' = 'treatment.ch1')
  
  # return as tibble to get rid of rownames
  return(tibble(md_filtered))

}

get_levels <- function(row){
  # helper function to return the levels of a row that is a factor
  row_levels <- levels(factor(row))
  return(paste(row_levels, collapse = ", "))
}

summarize_md <- function(md) {
  # create a summarized version of the metadata
  removes <- c('Sample name', 'GEO Accession')
    
  md_factors <- md %>%
    dplyr::select(!all_of(removes))  %>% # further restrict columns
    t() # transpose
  
  # create the summary
  md_summary <- tibble(
    'Column Name' = rownames(md_factors),
    'Type' = 'factor',
    'Distinct Values' = apply(md_factors, 1, get_levels)
  )
  
  return(md_summary)
  
}

visualize_md <- function(md, group, metric) {
  # visualize how the samples are broken down into different groups
  barplot <- ggplot(md, aes(fill=!!sym(metric), x=!!sym(group))) + 
    geom_bar(stat="count", 
             width=0.8) +
    labs(y=paste("Count of", metric),
         x=group,
         fill='',
         title=paste("Count of", metric, "per", group)) +
    theme_classic(base_size = 16) +
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(color='gray'),
          axis.line.x=element_line(color="gray"),
          axis.line.y=element_line(color="gray")) +
    scale_fill_brewer(palette="Dark2") +
    geom_text(aes(label = after_stat(count)), 
              stat = "count", 
              position = position_stack(vjust = 0.5),
              color="gray16")
  return(barplot)
}

get_data <- function(meta) {
  # divide samples up to prep for DESEQ2 analysis
  male_titles <- meta[meta$'Sex.ch1' == 'male','title'] 
  # remove the male outlier
  male_titles <- male_titles[! male_titles %in% c('CAL2M')]
  female_titles <- meta[meta$'Sex.ch1' == 'female','title']
  larvae_titles <- meta[meta$'lifestage.ch1' == 'larvae','title']
  
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

########################
######## DESEQ2 ######## 
########################

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
                                design = ~ timepoint + condition + timepoint:condition)
  
  # filter for rows with at least the number of reads indicated by count_filter
  keep <- rowSums(counts(dds)) >= count_filter
  dds <- dds[keep,]
  
  # do differential expression analysis
  dds <- DESeq(dds, test="LRT", reduced=~ timepoint + condition)
  
  # save results
  results <- results(dds, name=condition_name)

  return(results)
}

# get genes that are significant(padj < 0.01)
de_genes <- function(deseq_results){
  deseq_results <- deseq_results[!is.na(deseq_results$padj),]
  de_genes <- row.names(deseq_results[deseq_results$padj < 0.01,])
  return(de_genes)
}

