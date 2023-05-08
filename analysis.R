#!/usr/bin/Rscript
## Caroline Sullivan
## BF591 Final Project
## Analysis

libs <- c("gplots", "GEOquery", "tidyverse", "RColorBrewer", "ggVennDiagram", "BiocManager",
          "DESeq2", "DEGreport", "maSigPro")

for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}

################################# 
######## HELPER FUNCTIONS ####### 
################################# 

# CPM normalization
get_library_size <- function(count_data) {
  # transpose the tibble
  count_data_T <- count_data %>%
    pivot_longer(cols = -gene, names_to = 'sample') %>%
    pivot_wider(names_from = gene, values_from = value)
  # make new col with total counts
  count_data_T <- count_data_T %>%
    dplyr::mutate(library_size=rowSums(count_data_T[-1]))
  # return cols of interest
  return(dplyr::select(count_data_T, c(sample, library_size)))
}

normalize_by_cpm <- function(count_data) {
  # first get lib sizes and pivot
  lib_size <- get_library_size(count_data) %>% 
    pivot_wider(names_from=sample, values_from=library_size)
  
  # now use to normalize
  count_data_norm <- count_data
  count_data_norm[-1] <- mapply(function(x, y) x/y*1e6, x=count_data[-1], y=lib_size)
  
  return(count_data_norm)
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
  results <- list(male_counts, female_counts, larvae_counts, counts)
  names(results) <- c('Male', 'Female', 'Larvae', 'All')
  return(results)
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
  
  # normalize data with CPM method
  counts_data <- normalize_by_cpm(as_tibble(counts_data, rownames='gene'))
  counts_data <-as.data.frame(counts_data)
  rownames(counts_data) <- counts_data$gene
  counts_data <- subset(counts_data, select = -c(gene))
  counts_data <- prep_scatter_data(counts_data)
  
  # find the percentile cut-off for variance
  percentile_cutoff <- quantile(counts_data$variance, probs = slider_var)
  
  # create color vector
  color <- ifelse(counts_data['num_nonzero']>=slider_nonzero & counts_data['variance']>=percentile_cutoff, 
                  'TRUE', 
                  'FALSE')
  
  scatter1 <- ggplot(counts_data,
                     aes(x=log10(median_count + 1),y=log10(variance + 1),color=color)) +
    geom_point(size=1.5) +
    theme_bw(base_size = 16) +
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
  
  # normalize by CPM and prep data
  counts_data <- normalize_by_cpm(as_tibble(counts_data, rownames='gene'))
  counts_data <-as.data.frame(counts_data)
  rownames(counts_data) <- counts_data$gene
  counts_data <- subset(counts_data, select = -c(gene))
  counts_data <- prep_scatter_data(counts_data)
  
  # find the percentile cut-off for variance
  percentile_cutoff <- quantile(counts_data$variance, probs = slider_var)
  
  # create color vector
  color <- ifelse(counts_data['num_nonzero']>=slider_nonzero & counts_data['variance']>=percentile_cutoff, 
                  'TRUE', 
                  'FALSE')
  
  scatter2 <- ggplot(counts_data,
                     aes(x=log10(median_count + 1),y=num_zero,color=color)) +
    geom_point(size=1.5) +
    theme_bw(base_size = 16) +
    scale_color_manual(values = c('FALSE' = "lightblue", 'TRUE' = "darkblue")) +
    ggtitle('Relationship Between Number of Zero Samples and Median Count') +
    labs(color='Passes Filtering Conditions',
         x=expression("log"[10]*"Median Count"),
         y="Number of Zero Samples")
  
  return(scatter2)
  
}

plot_heatmap <- function(counts_data, meta, slider_var, slider_nonzero){
  
  # normalize counts data
  counts_data <- normalize_by_cpm(as_tibble(counts_data, rownames='gene'))
  counts_data <-as.data.frame(counts_data)
  rownames(counts_data) <- counts_data$gene
  counts_data <- subset(counts_data, select = -c(gene))
  
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
  counts_log_scale <- log2(counts_filtered + 1)
  
  # get colsidecolors
  meta <- meta %>%
    # dplyr::filter(title != 'CAL2M') %>% # remove outlier
    dplyr::mutate('color_bar' = case_when((Sex.ch1 == 'female') ~ 'orange',
                                          (Sex.ch1 == 'male') ~ 'purple',
                                          (Sex.ch1 == 'unknown') ~ 'gray'))
  
  # Add extra space to right of plot area; change clipping to figure
  par(mar=c(5, 4, 2, 6))
  
  # custom colors
  coul <- colorRampPalette(brewer.pal(11, "RdYlBu"))(25)
  
  gplots::heatmap.2(as.matrix(counts_log_scale), 
                    scale="row",
                    labCol=FALSE,
                    labRow=FALSE,
                    key=TRUE,
                    trace="none",
                    ColSideColors=meta$color_bar,
                    col=coul,
                    main='Log2-Scaled Counts')
  
  legend("bottomright", title = "Sample Group",
         legend=c("Male","Female", "Larvae"), 
         fill=c("purple","orange", "gray"),
         inset=c(-0.1,0.88),
         xpd=TRUE, cex=0.9,
         bty="n")

}

plot_pca <- function(counts_data, meta, dim1, dim2, color, shape){
  
  # get normalized version of counts matrix
  counts_data <- normalize_by_cpm(as_tibble(counts_data, rownames='gene'))
  counts_data <-as.data.frame(counts_data)
  rownames(counts_data) <- counts_data$gene
  counts_data <- subset(counts_data, select = -c(gene))
  
  # log scale the data
  counts_data <- log10(counts_data + 0.5)

  # perform PCA
  pca_results <- prcomp(t(counts_data), center=TRUE, scale=FALSE)
  
  # get variance explained
  # first get the variance
  variance <- pca_results$sdev**2
  
  # now turn into % explained variance
  percent_explained <- round((variance / sum(variance))*100)
  
  # pull out columns of interest and label rows
  md_cut <- meta %>% 
    dplyr::select("treatment.ch1", "Sex.ch1", "timepoint.ch1", "geo_accession", "title")
  
  # adjust row and column names
  colnames(md_cut) <- c('Treatment', 'Sex', 'Timepoint', "geo_accession", "title")
  
  row.names(md_cut) <- md_cut$title
  
  # ready to make plot
  pca_biplot <- pca_results$x %>%
    merge(md_cut, by = 'row.names') %>% # add metadata
    ggplot(aes(x=!!sym(paste('PC', dim1, sep='')), y=!!sym(paste('PC', dim2, sep='')), color=!!sym(color), shape=!!sym(shape))) +
    geom_point(size=2) +
    theme_classic(base_size = 14) +
    ggtitle(paste("PCA for", 'PC', dim1, 'vs.', 'PC', dim2, "CPM-normalized and log10-scaled counts")) +
    labs(x=paste('PC', dim1, ': ', percent_explained[dim1], '% variance', sep=""),
         y=paste('PC', dim2, ': ', percent_explained[dim2], '% variance', sep=""),
         color=color,
         shape=shape) +
    scale_shape_discrete(labels=c(unique(md_cut[shape]))) +
    scale_colour_brewer(palette="Dark2", labels=c(unique(md_cut[color])))
  
  return(pca_biplot)
  
}

########################
######## DESEQ2 ######## 
########################

# get the col_data to input into the DESeq2 analysis
get_col_data <- function(counts, md, dataset) {
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
  
  # define the order of the timepoints
  if (dataset == 'Larvae'){
    category_order <- c('high1', 'low1')
    
  } else {
    category_order <- c('high2', 'low2', 'high3', 'low3')
  }
  
  col_data <- col_data %>% arrange(factor(timepoint, levels = category_order))
  
  return(col_data[c('condition', 'timepoint')])
}

# run DESeq2
run_deseq <- function(count_dataframe, coldata, count_filter, condition_name) {
  # make deseqdataset object
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_dataframe[,rownames(coldata)]),
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

# get genes that are significant(padj < padj_cutoff)
de_genes <- function(deseq_results, padj_cutoff){
  deseq_results <- deseq_results[!is.na(deseq_results$padj),]
  de_genes <- row.names(deseq_results[deseq_results$padj < 1 * 10^padj_cutoff,])
  return(de_genes)
}

# create the volcano plot - similar to assignment 7
volcano <- function(dataf, x_name, y_name, slider, color1, color2){
  
  # make sure no na values for x or y
  dataf_filtered <- dataf %>% 
    drop_na(c(!!sym(x_name), !!sym(y_name))) %>%
    dplyr::filter(!!sym(y_name) > 0)
  
  # make color labels
  color <- ifelse(dataf_filtered[y_name] < 1 * 10^slider, 'TRUE', 'FALSE')
  
  volc <- dataf_filtered %>%
    ggplot(mapping=aes(x=!!sym(x_name), y=-log10(!!sym(y_name)), color=color)) +
    geom_point(size=1.5) +
    theme_bw() + 
    scale_color_manual(values = c('FALSE' = color1, 'TRUE' = color2)) +
    labs(color=paste(y_name, '< 1e', slider),
         x=x_name,
         y=paste("-log(", y_name, ")", sep='')) +
    theme(legend.position="bottom")
  
  return(volc)
}

########################
###### CLUSTERING ######
########################

# perform clustering analysis on counts data using degpatterns
clustering <- function(dataset, counts, meta, deseq_results, padj_cutoff, minc){
  # get col data
  col_data <- get_col_data(counts, meta, dataset)
  
  # get the significant genes from the LRT deseq results at the indicated p value
  sigLRT_genes <- de_genes(deseq_results, padj_cutoff)
  
  # need to get rlog version of counts matrix
  rlog_mat <- DESeq2::vst(as.matrix(counts))
  
  # get the rlog counts value for the significant genes
  cluster_rlog <- rlog_mat[rownames(rlog_mat) %in% sigLRT_genes, ]
  
  # re-order the columns for graphing
  cluster_rlog <- cluster_rlog[,rownames(col_data)]
  
  # get clusters
  clusters <- degPatterns(cluster_rlog, 
                          metadata = col_data, 
                          time = "timepoint", 
                          col="condition",
                          scale=TRUE,
                          minc=minc)
  
  return(clusters)
}

# make plots for each of the clusters tracking gene expression over time
cluster_full_plot <- function(clusters){
  # set degree and level according to dataset
  if (length(levels(clusters$plot$data$timepoint)) == 2){
    degree <- 1
    level <- c('high1', 'low1')
  } else {
    degree <- 3
    level <- c('high2', 'low2', 'high3', 'low3')
  }
  
  
  # get the plotting data
  plot_data <- clusters$plot$data
  
  plot <- plot_data %>%
    ggplot(mapping = aes(x=factor(timepoint,level=level), 
                         y=value, 
                         color=condition,
                         fill=condition)) +
    
    geom_boxplot(outlier.alpha = NULL, alpha=0.3, outlier.size=0.9) +
    
    geom_point(na.rm=FALSE, position=position_jitterdodge(), alpha=0.3, size=0.9) +
    
    geom_smooth(mapping=(aes(x=factor(timepoint,level=level), 
                             y=value, 
                             color=condition,
                             group=condition)),
                method = lm, formula = y ~ poly(x, degree),
                se = FALSE,
                alpha=0.5) +
    
    geom_line(mapping=aes(group=line_group), alpha=0.1) +
    
    labs(x="",
         y='Z-score of gene abundance') +

    theme_bw() +
    scale_colour_brewer(palette="Set1")
    
  
  return(plot + facet_wrap('~title'))
}

# get the relavent data for download
get_summary_data <- function(clusters){
  # get the summary data
  summary <- clusters$normalized
  
  # keep relevant columns
  keep <- c('genes', 'merge', 'value', 'condition', 'timepoint', 'cluster')
  summary <- summary[,keep]
  
  return(summary)
  
}

# make the timecourse heatmap
timecourse_heatmap <- function(clusters){
  
  # get metadata and z score gene abundance values
  md <- clusters$normalized[,c('genes', 'merge', 'value', 'condition', 'timepoint', 'cluster')]
  plot_data <- md[,c('genes','merge', 'value')]
  plot_data <- plot_data %>% pivot_wider(names_from='merge', values_from='value')

  # make genes rownames
  plot_data <- as.data.frame(plot_data)
  rownames(plot_data) <- plot_data$genes
  plot_data <- subset(plot_data, select=-c(genes))

  # re-order the columns for graphing and get colsidecolors
  if (length(levels(clusters$plot$data$timepoint)) == 2){
    order <- c('controlhigh1', 'controllow1', 
               'fluctuatinghigh1', 'fluctuatinglow1')
    col_colors <- rep(c('green','red'),each=2)
    
  }else{
    order <- c('controlhigh2', 'controllow2', 'controlhigh3', 
               'controllow3','fluctuatinghigh2', 'fluctuatinglow2', 
               'fluctuatinghigh3', 'fluctuatinglow3')
    col_colors <- rep(c('green','red'),each=4)
  }
  
  # get row order
  row_order_df <- arrange(clusters$df, cluster)
  
  # get rows and columns in order
  plot_data <- plot_data[,order]
  plot_data <- plot_data[row_order_df$genes,]
  
  # get row colors
  coul <- brewer.pal(12, "Paired") 
  coul <- colorRampPalette(coul)(max(row_order_df$cluster))
  # randomize
  coul <- sample(coul)
  row_order_df$colors <- coul[row_order_df$cluster]
  
  # Add extra space to right of plot area; change clipping to figure
  par(mar=c(5.1, 4.2, 4.1, 8.4), xpd=TRUE)
  
  # have timepoints show up on heatmap
  colnames(plot_data) <- order
  
  gplots::heatmap.2(as.matrix(plot_data), 
                    Colv=NULL,
                    Rowv=NULL,
                    dendrogram="none",
                    labRow=FALSE,
                    key=TRUE,
                    trace="none",
                    ColSideColors=col_colors,
                    RowSideColors=row_order_df$colors,
                    srtCol=28,
                    col=brewer.pal(11, "RdYlBu"),
                    main='z-score of gene abundance by cluster')
  
  legend("bottomright", title = "Condition",
         legend=c("control","fluctuating"), 
         fill=c("green","red"),
         inset=c(-0.1,0.88),
         xpd=TRUE, cex=0.9,
         bty="n")
  
  legend("bottomright", title="Clusters",
         legend=unique(row_order_df[c('cluster', 'colors')])$cluster,
         fill=unique(row_order_df[c('cluster', 'colors')])$colors,
         inset=c(0.88,0.1),
         xpd=TRUE, cex=0.9,
         bty="n")
  
}

# show z-score gene abundance for individual gene of interest
# make plots for each of the clusters tracking gene expression over time
cluster_gene_plot_agg <- function(clusters, dataset, gene_name){
  # set degree and level according to dataset
  if (dataset == "Larvae"){
    level <- c('high1', 'low1')
  } else {
    level <- c('high2', 'low2', 'high3', 'low3')
  }
  
  
  # get the plotting data and subset for the specific gene of interest
  plot_data <- clusters$plot$data
  plot_data <- plot_data %>%
    dplyr::filter(genes==gene_name)
  
  plot <- plot_data %>%
    ggplot(mapping = aes(x=factor(timepoint,level=level), 
                         y=value, 
                         color=condition,
                         fill=condition)) +
    
    geom_point(na.rm=FALSE, alpha=0.7, size=0.9) +
    
    geom_line(mapping=aes(group=condition), alpha=0.7) +
    
    labs(x="",
         y='Z-score of gene abundance',
         title=paste('Expression Over Time for', gene_name)) +
    
    theme_bw() +
    scale_colour_brewer(palette="Set1")
  
  
  return(plot)
}


# alt version - uses individual counts from each replicate
cluster_gene_plot <- function(clusters, col_dat, gene_name){
  # set degree and level according to dataset
  if (length(levels(clusters$plot$data$timepoint)) == 2){
    level <- c('high1', 'low1')
  } else {
    level <- c('high2', 'low2', 'high3', 'low3')
  }
  
  
  # get the plotting data and subset for the specific gene of interest
  plot_data <- clusters$counts
  plot_data <- plot_data[gene_name,]
  
  # add the metadata
  plot_data <- cbind(plot_data, col_dat)
  colnames(plot_data) <- c('value', 'condition', 'timepoint')
  plot_data$line_group <- paste(plot_data$condition, plot_data$timepoint, sep="_")
  
  plot <- plot_data %>%
    ggplot(mapping = aes(x=factor(timepoint,level=level), 
                         y=value, 
                         color=condition,
                         fill=condition)) +
    
    geom_boxplot(outlier.alpha = NULL, alpha=0.5, outlier.size=0.9) +
    
    geom_point(na.rm=FALSE, position=position_jitterdodge(), alpha=0.5, size=0.9) +
    
    geom_line(mapping=aes(group=condition), alpha=0.7) +
    
    labs(x="",
         y='vst normalized gene expression values',
         title=paste('Expression Over Time for', gene_name)) +
    
    theme_bw() +
    scale_colour_brewer(palette="Set1")
  
  
  return(plot)
}

# ALTERNATIVE OPTION: MASIGPRO
# get the right design format to run masigpro by editing coldata
make_masig_design <- function(col_dat){
  # make separate columns for timepoints
  col_dat$Time <- as.numeric(col_dat$timepoint)
  
  # separate columns for replicates
  col_dat <- col_dat %>% dplyr::mutate(reps = paste(condition, timepoint))
  col_dat$Replicate <- as.numeric(as.factor(col_dat$reps))
  
  # separate columns for conditions
  col_dat <- col_dat %>%
    dplyr::mutate('Control' = ifelse(condition == 'control', 1, 0), 
                   'Fluctuating' = ifelse(condition == 'fluctuating', 1, 0))

  return(col_dat[c('Time', 'Replicate', 'Control', 'Fluctuating')])
  
}
  
# perform clustering analysis on counts data using masigpro as seen in the paper
clustering_masigpro <- function(counts, coldat, degree, k){
  
  # get design data
  des <- make_masig_design(coldat)
  
  # make design
  design <- make.design.matrix(des, degree = degree)
  
  # need to get normalized version of counts matrix
  counts_norm <- normalize_by_cpm(as_tibble(counts, rownames='gene'))
  counts_norm <-as.data.frame(counts_norm)
  rownames(counts_norm) <- counts_norm$gene
  counts_norm <- subset(counts_norm, select = -c(gene))
  
  # get significant genes
  fit <- p.vector(counts_norm, 
                  design, 
                  Q = 0.05, 
                  MT.adjust = "BH", 
                  min.obs = 10, 
                  counts=TRUE)
  
  # do stepwise regression
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
  
  # get signficant genes
  sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
  
  par(mar = c(1, 1, 1, 1))
  
  see.genes(sigs$sig.genes$FluctuatingvsControl, show.fit = T, dis =design$dis,
            cluster.method="kmeans", cluster.data = 1, k = k)
  
}

