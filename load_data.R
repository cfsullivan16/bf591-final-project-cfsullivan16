#!/usr/bin/Rscript
## Caroline Sullivan
## BF591 Final Project
## Loading in the data

libs <- c("GEOquery", "tidyverse")

for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}

# GET DATA
#gse <- GEOquery::getGEO("GSE64810", GSEMatrix = TRUE)
#meta <- pData(gse[[1]])
#filePaths = GEOquery::getGEOSuppFiles("GSE64810")

# unzip the files
#GEOquery::gunzip('~/Documents/BF591/bf591-final-project-cfsullivan16/GSE64810/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz', 
#                 destname = gsub("[.]gz$", "", '~/Documents/BF591/bf591-final-project-cfsullivan16/data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz'), 
#                 overwrite = FALSE,
#                 remove = TRUE, 
#                 BFR.SIZE = 1e+07)

#GEOquery::gunzip('~/Documents/BF591/bf591-final-project-cfsullivan16/GSE64810/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz', 
#                 destname = gsub("[.]gz$", "", '~/Documents/BF591/bf591-final-project-cfsullivan16/data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz'), 
#                 overwrite = FALSE,
#                 remove = TRUE, 
#                 BFR.SIZE = 1e+07)

#filePaths2 = GEOquery::getGEOSuppFiles("GSE64403")

#GEOquery::gunzip('~/Documents/BF591/bf591-final-project-cfsullivan16/GSE64403/GSE64403_FPKM_table.txt.gz', 
#                 destname = gsub("[.]gz$", "", '~/Documents/BF591/bf591-final-project-cfsullivan16/data/GSE64403_FPKM_table.txt.gz'), 
#                 overwrite = FALSE,
#                 remove = TRUE, 
#                 BFR.SIZE = 1e+07)

# get the counts matrix and unzip
filePaths3 = GEOquery::getGEOSuppFiles("GSE150450")
GEOquery::gunzip('~/Documents/BF591/bf591-final-project-cfsullivan16/GSE150450/GSE150450_gene_count_matrix.csv.gz', 
                 destname = gsub("[.]gz$", "", '~/Documents/BF591/bf591-final-project-cfsullivan16/data/GSE150450_gene_count_matrix.csv.gz'), 
                 overwrite = FALSE,
                 remove = TRUE, 
                 BFR.SIZE = 1e+07)

# get the metadata
gse3 <- GEOquery::getGEO("GSE150450", GSEMatrix = TRUE)
meta3 <- pData(gse3[[1]])
write.csv(meta3, 'data/meta.csv', row.names=TRUE)

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

male_counts <- counts[,male_titles$title]
female_counts <- counts[,female_titles]
larvae_counts <- counts[,larvae_titles$title]

##### DESEQ2 for FEMALE SAMPLES #####


