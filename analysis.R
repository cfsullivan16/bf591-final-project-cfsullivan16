#!/usr/bin/Rscript
## Caroline Sullivan
## BF591 Final Project
## Analysis

libs <- c("gplots", "GEOquery", "tidyverse", "RColorBrewer", "ggVennDiagram", "BiocManager",
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

################################
### FUNCTIONS FOR COUNTS TAB ###
################################

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
  counts <- subset(counts, select = -c(gene_id))
  
  male_counts <- counts[,male_titles]
  female_counts <- counts[,female_titles]
  larvae_counts <- counts[,larvae_titles]
  
  # put all counts back together, excluding male outlier
  # all_counts <- cbind(male_counts, female_counts, larvae_counts)
  
  return(list(male_counts, female_counts, larvae_counts, counts))
}

filter_counts <- function(counts_data, slider_var, slider_nonzero){
  # takes in filtering conditions and spits out a dataframe containing
  # summary info about filtering results
  
  # get initial number of samples and genes
  num_samps <- ncol(counts_data)
  num_genes <- nrow(counts_data)
  
  # first, calculate variance and number of non-zero samples
  variance <- apply(counts_data, 1, var)
  num_zero <- rowSums(counts_data!=0)
  
  counts_data <- counts_data %>% 
    dplyr::mutate('variance'=variance,
                  'num_nonzero'=num_zero)
  
  # find the percentile cut-off for variance
  percentile_cutoff <- quantile(counts_data$variance, probs = slider_var)
  
  # now apply the filters
  counts_filtered <- counts_data %>%
    dplyr::filter(variance >= percentile_cutoff,
                  num_nonzero >= slider_nonzero) %>%
    dplyr::select(-c('variance', 'num_nonzero')) # can get rid of these cols now
  
  # calculate number that passed and failed
  num_genes_passed <- nrow(counts_filtered)
  num_genes_failed <- num_genes - num_genes_passed
  
  # make the final table
  sum_table <- tibble(
    "Number of Samples" = scales::comma(num_samps),
    "Number of Genes" = scales::comma(num_genes),
    "Genes Passing Current Filter" = paste(scales::comma(num_genes_passed), 
                                           " (",
                                           round((num_genes_passed/num_genes)*100, digits=2),
                                           "%)", sep=""),
    "Genes NOT Passing Current Filter" = paste(scales::comma(num_genes_failed),
                                               " (",
                                               round((num_genes_failed/num_genes)*100, digits=2),
                                               "%)", sep="")
  )
  
  # transpose the table for readability
  sum_table <- t(as.data.frame(sum_table))
  colnames(sum_table) = 'Result'
  sum_table <- sum_table %>% as_tibble(rownames='Metric')
  
  return(sum_table)
  
}

prep_scatter_data <- function(counts_data){
  # get necessary columns to make both scatter plots
  # first, calculate variance and number of zero samples and nonzero samples
  variance <- apply(counts_data, 1, var)
  num_zero <- rowSums(counts_data==0)
  num_nonzero <- rowSums(counts_data!=0)
  
  # get median, then add other cols
  counts_data <- counts_data %>% 
    dplyr::mutate('median_count'=apply(., 1, median, na.rm=TRUE)) %>%
    dplyr::mutate('variance'=variance,
                  'num_zero'=num_zero,
                  'num_nonzero'=num_nonzero)
  
  return(counts_data)
  
}

scatter_plot1 <- function(counts_data, slider_var, slider_nonzero){
  # diagnostic scatter plots - genes passing filters marked in a darker color,
  # genes filtered out are lighter
  # scatter_plot1: median count vs variance (consider log scale for plot)
  # scatter_plot2: median count vs number of zeros
  counts_data <- prep_scatter_data(counts_data)
  
  # find the percentile cut-off for variance
  percentile_cutoff <- quantile(counts_data$variance, probs = slider_var)
  
  # create color vector
  color <- ifelse(counts_data['num_nonzero']>=slider_nonzero & counts_data['variance']>=percentile_cutoff, 
                  'TRUE', 
                  'FALSE')
  
  scatter1 <- ggplot(counts_data,
                     aes(x=log10(median_count),y=log10(variance),color=color)) +
    geom_point(size=1.5) +
    theme_bw() +
    scale_color_manual(values = c('FALSE' = "lightblue", 'TRUE' = "darkblue")) +
    ggtitle('Relationship Between Variance and Median Count') +
    labs(color='Passes Filtering Conditions',
         x=expression("log"[10]*"Median Count"),
         y=expression("log"[10]*"Variance"))
  
  return(scatter1)
  
}

scatter_plot2 <- function(counts_data, slider_var, slider_nonzero){
  # diagnostic scatter plots - genes passing filters marked in a darker color,
  # genes filtered out are lighter
  # scatter_plot1: median count vs variance (consider log scale for plot)
  # scatter_plot2: median count vs number of zeros
  counts_data <- prep_scatter_data(counts_data)
  
  # find the percentile cut-off for variance
  percentile_cutoff <- quantile(counts_data$variance, probs = slider_var)
  
  # create color vector
  color <- ifelse(counts_data['num_nonzero']>=slider_nonzero & counts_data['variance']>=percentile_cutoff, 
                  'TRUE', 
                  'FALSE')
  
  scatter2 <- ggplot(counts_data,
                     aes(x=log10(median_count),y=num_zero,color=color)) +
    geom_point(size=1.5) +
    theme_bw() +
    scale_color_manual(values = c('FALSE' = "lightblue", 'TRUE' = "darkblue")) +
    ggtitle('Relationship Between Variance and Number of Zero Samples') +
    labs(color='Passes Filtering Conditions',
         x=expression("log"[10]*"Median Count"),
         y="Number of Zero Samples")
  
  return(scatter2)
  
}

plot_heatmap <- function(counts_data, meta, slider_var, slider_nonzero){
  # get information we need to filter out genes
  counts_data <- prep_scatter_data(counts_data)
  
  # find the percentile cut-off for variance
  percentile_cutoff <- quantile(counts_data$variance, probs = slider_var)
  
  # now apply the filters
  counts_filtered <- counts_data %>%
    dplyr::filter(variance >= percentile_cutoff,
                  num_nonzero >= slider_nonzero) %>%
    dplyr::select(-c('variance', 'num_nonzero', 'num_zero', 'median_count')) # can get rid of these cols now
  
  # log scale the counts
  counts_log_scale <- log10(counts_filtered + 1)
  
  # get rowsidecolors
  meta <- meta %>%
    # dplyr::filter(title != 'CAL2M') %>% # remove outlier
    dplyr::mutate('color_bar' = case_when((Sex.ch1 == 'female') ~ 'orange',
                                          (Sex.ch1 == 'male') ~ 'purple',
                                          (Sex.ch1 == 'unknown') ~ 'gray'))
  
  
  return(gplots::heatmap.2(as.matrix(counts_log_scale), 
                           scale="row",
                           labCol=FALSE,
                           labRow=FALSE,
                           key=TRUE,
                           trace="none",
                           ColSideColors=meta$color_bar,
                           col=brewer.pal(11, "RdYlBu")))
}

plot_pca <- function(counts_data, meta, slider_var, slider_nonzero){
  # log scale the data
  counts_data <- log10(counts_data + 1)

  # perform PCA
  pca_results <- prcomp(t(counts_data), center=TRUE, scale=FALSE)
  
  # pull out columns of interest and label rows
  md_cut <- meta %>% 
    dplyr::select("treatment.ch1", "Sex.ch1", "geo_accession", "title")
  
  row.names(md_cut) <- md_cut$title
  
  # ready to make plot
  pca_biplot <- pca_results$x %>%
    merge(md_cut, by = 'row.names') %>% # add metadata
    ggplot(aes(x=PC1, y=PC2, color=treatment.ch1, shape=Sex.ch1)) +
    geom_point() +
    scale_colour_brewer(palette="Dark2") +
    theme_classic() 
  
  return(pca_biplot)
  
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

