#library('sys')
library('DESeq2')
library('ggplot2')
#library('stringr')
#library('GenomicRanges')
library('dplyr')

source('/home/cwake/snakemakes/Utility_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')
print(sessionInfo())

if(interactive()){
  project <- '2021612_finch'
  test <- 'Treatment'
  adjust <- 'Sex'
  batch <- 'plate, RNA_Extraction_Batch'
  qc_name <- '2024-03-08'
  # project <- '2022612_Petrovas'
  # test <- 'Strain'
  # adjust <- 'Cell_type'
  # batch <- ''
  # qc_name <- 'QC3'
  # 
  count_file <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/filteredCounts.txt')
  covs_file <- paste0('/home/cwake/projects/', project, '/data/Covariates_QC_metrics_intermediate.csv')
  #covs_file <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_filter.csv')
  norm_file  <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/Treatment/normalizedCounts.txt')
  dds_file  <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/Treatment/normalizedCounts.RDS')
  pdf_out <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/normalization.pdf')
} else{
  args = commandArgs(trailingOnly=TRUE)
  count_file <- args[1]
  covs_file <- args[2]
  test <- args[3]
  adjust <- args[4]
  batch <- args[5]
  norm_file <- args[6]
  dds_file <- args[7]
  pdf_out <- args[8]
}

covs_names <- c(test, 'ID', 'Sex')
names(covs_names) <- c('Test', 'ID', 'Sex')

test <- trimws(strsplit(test,',')[[1]])[1]
adjust <- trimws(strsplit(adjust,',')[[1]])
batch <- trimws(strsplit(batch,',')[[1]])
### Read counts
raw_counts <- read.csv(count_file, stringsAsFactors=F, header = T,
                       check.names = F, 
                       colClasses = 'character',
                       sep = '\t')
### Convert to numeric
raw_counts <- mutate_all(raw_counts, function(x) as.numeric(x))

### Read in covariates
covs <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
### Remove the filtered samples
covs <- covs[which(covs$filter == '0'),]
row.names(covs) <- covs$Sample_ID
raw_counts <- raw_counts[, row.names(covs)]

all_covs <- c('ID', test, adjust, batch)
all_covs <- all_covs[which(all_covs != '')]
### model_covs_deseq is an input array
deseq_model <- Make_DESeq2_model(test, adjust)

### Convert covs data frame to DataFrame (S4) object, as necessary for DESeq2
covs <- DataFrame(covs[match(colnames(raw_counts), covs[,'ID']) ,])
covs <- covs[, all_covs]
for(col in colnames(covs)){
  covs[,col] <- gsub('-', '', covs[,col])
  ### This is really just a best guess at data type based on name of the column, uniqueness, and numeric conversion
  covs[,col] <- determine_type(covs, col)
}

### DESeq2 Normalization, returning norm_counts and dds (DESeqDataSet)
output <- DE_normalization(test, deseq_model, raw_counts, covs, do_transform = FALSE)
dds <- output[[1]]
norm_counts <- data.frame(output[[2]], check.names = F)


plotdat <- as.data.frame(colSums(raw_counts))
plotdat$data <- 'raw'
colnames(plotdat) <- c('Nreads', 'transformation')
plotdat2 <- as.data.frame(colSums(norm_counts))
plotdat2$data <- 'normalized'
colnames(plotdat2) <- c('Nreads', 'transformation')
plotdat <- rbind(plotdat, plotdat2)
p <- ggplot(plotdat, aes(x = Nreads, color = transformation)) +
  geom_histogram(fill = "white", position = "dodge")+
  theme(legend.position = "top")
pdf(pdf_out)
p
dev.off()

write.table(norm_counts, norm_file, quote = F, sep = '\t', row.names = T, col.names = T)
saveRDS(dds, file = dds_file)
