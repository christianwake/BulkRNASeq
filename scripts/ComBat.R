#library('sys')
#library('DESeq2')
#library('ggplot2')
#library('stringr')
#library('GenomicRanges')
#library('dplyr')
library('readr')

source('/home/cwake/snakemakes/Utility_functions.R')
#source('/home/cwake/snakemakes/DE_functions.R')

if(interactive()){
  #project <- '2021612_finch'
  project <- '2022612_Petrovas'
  count_in <- paste0('/home/cwake/projects/', project, '/data/counts/normalizedCounts.txt')
  dds_in <- paste0('/home/cwake/projects/', project, '/data/counts/normalizedCounts.RDS')
  covs_in <- paste0('/home/cwake/projects/', project, '/data/Covariates_QC_metrics_filter.csv')
  batch_info <- paste0('/home/cwake/projects/', project, '/data/batch_evaluation.txt')
  count_out <- paste0('/home/cwake/projects/', project, '/data/counts/finalCounts.txt')
  dds_out <- paste0('/home/cwake/projects/', project, '/data/counts/finalCounts.RDS')
} else{
  args = commandArgs(trailingOnly=TRUE)
  count_in <- args[1]
  dds_in <- args[2]
  covs_in <- args[3]
  batch_info <- args[4]
  count_out <- args[5]
  dds_out <- args[6]
}

### Read counts
norm_counts <- read.csv(count_in, stringsAsFactors=F, header = T,
                        check.names = F, 
                        colClasses = 'character',
                        sep = '\t')
### Convert to numeric
norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))
dds <- readRDS(dds_in)

### Read in covariates
covs <- read.csv(covs_in, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
covs <- covs[which(covs$filter == '0'),]

batch <- read_file(batch_info)
if(grepl('Enter batch information', batch) | (!grepl('Combat batch', batch))){
  print('No ComBat batch correction is performed.')
  write.table(norm_counts, count_out, quote = F, sep = '\t', row.names = T, col.names = T)
  saveRDS(dds, dds_out)
} else{
  batch <- strsplit(batch, '\n')[[1]]
  combat <- batch[[which(sapply(1:length(batch), function(x) grepl('ComBat batch', batch[x])))]]
  combat <- gsub('ComBat batch: ', '', combat)
  print('Placeholder.')
  write.table(norm_counts, count_out, quote = F, sep = '\t', row.names = T, col.names = T)
  saveRDS(dds, dds_out)
}

