#library('sys')
library('DESeq2')
library('ggplot2')
#library('stringr')
#library('GenomicRanges')
library('dplyr')
library('readr')
library('effsize')
library('data.table')
library('biomaRt')
library('EnhancedVolcano')

source('/home/cwake/snakemakes/Utility_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')

print(Sys.Date())

if(interactive()){
  project <- '2021612_finch'
  #qc_name <- '2023-12-06'
  qc_name <- '2024-03-08'
  test <- 'Treatment'
  adjust <- 'Sex,Euth_Age,plate'
  gtf_file <- '/home/cwake/resources/genomes/tguttata/bTaeGut1_v1.p/Annotation/Taeniopygia_guttata.bTaeGut1_v1.p.104.gtf'
  stratification <- 'All'

  count_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/finalCounts.txt')
  dds_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/finalCounts.RDS')
  covs_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_filter.csv')
  batch_info <- paste0('/home/cwake/projects/', project, '/data/batch_evaluation.txt')
  pdf_out <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test,  '/', stratification, '/DESeq2_results.pdf')
  txt_out <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/', stratification, '/DESeq2_results.txt')
}else{
  args = commandArgs(trailingOnly=TRUE)
  count_in <- args[1]
  dds_in <- args[2]
  covs_in <- args[3]
  batch_info <- args[4]
  gtf_file <- args[5]
  test <- args[6]
  adjust <- args[7]
  stratification <- args[8]
  txt_out <- args[9]
  pdf_out <- args[10]  
}

adjust <- trimws(strsplit(adjust, ',')[[1]])

### Read batch file and add to model covariates
#scan(batch_info, what = character(), quote = "", sep = ':', strip.white = T)
batch <- read_file(batch_info)
batch <- strsplit(batch, '\n')[[1]]
batch <- batch[[which(sapply(batch, function(x) grepl('Model batch', x)))]]
batch <- gsub('Model batch: ', '', batch)
model_covs <- c(adjust, batch)
### In case batch is empty
model_covs <- model_covs[which(model_covs != '')]

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
row.names(covs) <- covs$ID
### Filter samples from covs and counts, and the equivalent in dds
covs <- covs[which(covs$filter == '0'),]
norm_counts <- norm_counts[, colnames(norm_counts)[which(colnames(norm_counts) %in% row.names(covs))]]
dds <- dds[, covs$ID]

colkeep <- c('ID', test, model_covs)
if(!is.na(stratification) & stratification != '' & stratification != 'All'){
  colkeep <- c(colkeep, strsplit(stratification, '-')[[1]][1])
}

### Convert covs data frame to DataFrame (S4) object, as necessary for DESeq2
covs <- DataFrame(covs)
covs <- covs[, colkeep]
for(col in colnames(covs)){
  #covs[, col] <- droplevels(covs[, cov])
  ### This is really just a best guess at data type based on name of the column, uniqueness, and numeric conversion
  covs[, col] <- determine_type(covs, col)
}
### In case an entire factor level is now unused after filtering, droplevels()
for(cov in model_covs){
  #covs[, cov] <- droplevels(covs[, cov])
  dds@colData[, cov] <- droplevels(dds@colData[, cov])
}
##### Parse stratification input
### If stratification is input (not 'All', NA or '')
if(stratification != '' & !is.na(stratification) & stratification != 'All'){
  a <- strsplit(stratification, '-')[[1]]
  strat <- a[2]
  names(strat) <- a[1]
  ### If the input stratification actually matches a column name
  if(names(strat) %in% colnames(covs)){
    ### If this subset has fewer than 5 samples
    if(table(covs[, names(strat)])[strat] < 5){
      stop('Cannot do stratification. Very few samples in subset.')
    }
    ### Remove the stratification cov from adjust, if there.
    adjust <- adjust[which(!adjust %in% names(strat))]
    ### Remove samples from covs, norm_counts, and dds
    covs <- covs[which(covs[, names(strat)] == strat), ]
    ### Option 1 - remove from norm_counts and dds, and model explicicately changed with design() assignment after batch chunk
    norm_counts <- norm_counts[, covs$ID]
    dds <- dds[, covs$ID]
    ### In case an entire factor level is now unused, droplevels()
    for(cov in model_covs){
      covs[, cov] <- droplevels(covs[, cov])
      dds@colData[, cov] <- droplevels(dds@colData[, cov])
    }
    ### Option 2 - Renormalize???
    #raw_counts <- raw_counts[, covs$ID]
    # deseq_model <- Make_DESeq2_model(test, adjust)
    # output <- DE_normalization(test, deseq_model, raw_counts, covs, do_transform = FALSE)
    # dds <- output[[1]]
    # norm_counts <- data.frame(output[[2]], check.names = F)
    
  } else{
      stop('Cannot do stratification. Input does not match a covariate name')
  }
}

model_covs <- c(adjust, batch)
model_covs <- model_covs[which(model_covs != '')]

###### Set up the DESeq2 and Firth model formats
### Create the correct DESeq2 model format as a string
deseq_model <- Make_DESeq2_model(test, model_covs)
### Manually change model?
design(dds) <- eval(parse(text=deseq_model))

print('Model:')
print(deseq_model)
###################
###### Run DESeq2
###################
print("DESeq2")
### Outlier trimming  by default identifies outlier counts of genes with more than 7 sample replicates, and replaces outliers by the mean. This is difficult to do if you have a quantitative covariate in the model.
### Independent filtering removes features from the testing that have little to no chance to show a positive result.
deseq_output <- Run_DESeq2(dds, test, model_covs, covs, pdf_file = pdf_out, txt_file = NA, outlier_trim = F, independent_filtering = T, cohen = F)

# norm_counts <- data.frame(deseq_output[1], check.names = F)
# write.table(norm_counts, file = paste0(DE_path,"/DESeq2_normalized_counts.txt"), quote = F, sep = ',')
res <- deseq_output[2][[1]]
res_df <- as.data.frame(res)

print('Reading gtf file')
print(gtf_file)
if(gtf_file == '' | is.na(gtf_file)){
  gtf <- NA
} else{
  if(grepl('\\.gtf', gtf_file)){
    gtf <- read_gtf(gtf_file, 'gene', c('gene_name', 'gene_id')) 
  } else{
    gtf <- readRDS(gtf_file)
  }
  row.names(gtf) <- gtf$gene_id
  ### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
  gtf_cols <- c('gene_name', 'gene_id')
  #id_type <- names(which.max(sapply(gtf_cols, function(col) sum(row.names(norm_counts) %in% gtf[, col]))))
  id_type <- get_id_name(row.names(norm_counts), gtf)
}
res[, 'gene_name_gtf'] <- gtf[row.names(res), 'gene_name']

### Gene name first for legibility
res <- res[, c('gene_name_gtf', head(colnames(res), n = length(colnames(res))-1))]

### Order by p-value
res <- res[order(res$pvalue), ]
write.table(res, file = txt_out, quote = F, sep = '\t')
