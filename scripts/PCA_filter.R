#library('sys')
library('DESeq2')
library('ggplot2')
library('stringr')
library('GenomicRanges')
library('dplyr')
library('ggrepel')

source('/home/cwake/snakemakes/Utility_functions.R')
#source('/home/cwake/snakemakes/DE_functions.R')

if(interactive()){
  project <- '2021612_finch'
  #project <- '2022612_Petrovas'
  qc_name <- '2023-12-06'
  count_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/normalizedCounts.txt')
  covs_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_intermediate.csv')
  qc_file <- paste0('/home/cwake/projects/', project, '/QC_steps/Sample_and_feature_filters.csv')
  covs_out <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_filter.csv')
  out_pdf <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/PCA_filter.pdf')
  out_csv <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_filter.csv')
  
} else{
  args = commandArgs(trailingOnly=TRUE)
  count_in <- args[1]
  covs_in <- args[2]
  qc_file <- args[3]
  covs_out <- args[4]
  out_pdf <- args[5]
  out_csv <- args[6]
}

filters <- read.table(qc_file, header = T, sep = ',')

### ignore everything except PCA
filters <- filters[which(filters$feature == 'PCA'),]

pca <- filters[which(filters$type  == 'sample'), 'value']

### Read counts
norm_counts <- read.csv(count_in, stringsAsFactors=F, header = T,
                       check.names = F, 
                       colClasses = 'character',
                       sep = '\t')
### Convert to numeric
norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))

### Read in covariates
covs <- read.csv(covs_in, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')

### Set default
if(is.na(pca) | pca == ''){
  pca <- '3(SD)'
}
pca_sd <- as.numeric(gsub('\\(SD\\)', '', pca)) ### PCA thresholds can only be in units of standard deviations from the mean

PCA <- run_pca(norm_counts)
### Returned outlier values are the fraction of a SD above/below the mean
outliers <- pca_outliers(PCA, c(1,2,3), pca_sd)

### Update covs
covs[which(covs$ID %in% row.names(outliers[which(rowSums(outliers != 'Not') > 0),])), 'filter'] <- 'PCA'
write.table(covs, covs_out, quote = F, sep = ',', row.names = F, col.names = T)

outs <- outliers[which(rowSums(outliers != 'Not') == 1), ]
write.table(outs, out_csv, quote = F, sep = ',', row.names = T, col.names = T)

outliers[outliers != 'Not'] <- 'Outlier'
pdf(out_pdf)
pca_plot_outliers(PCA, outliers, c(1,2), note = paste0(', outliers ', pca))
pca_plot_outliers(PCA, outliers, c(1,3), note = paste0(', outliers ', pca))
dev.off()



