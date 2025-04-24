library('sys')
library('ggplot2')
library('dplyr')

source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021612_finch'
  project <- '2021612_finch/RRBS'
  qc_name <- '2023-04-13'
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/All_covariates.csv')
  count_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/counts/featureCounts.txt')
  star_path <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/bam/')
  out_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/Covariates_QC_metrics.csv')
} else{
  args = commandArgs(trailingOnly=TRUE)
  
  covs_file <- args[1]
  star_path <- args[2]
  count_file <- args[3]
  out_file <- args[4]
  
}

### Read in covariates
covs <- read.csv(covs_file,stringsAsFactors=F,header = T,
                 check.names=F, 
                 colClasses = 'character',
                 sep=',')

### Read results of featureCounts
raw_counts <- read.csv(count_file, stringsAsFactors=F, header = T,
                       check.names=F, 
                       colClasses = 'character',
                       sep='\t', skip = 1)
row.names(raw_counts) <- raw_counts[, 'Geneid']
raw_counts <- raw_counts[, which(! colnames(raw_counts) %in% c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'))]
### Convert to numeric
raw_counts <- mutate_all(raw_counts, function(x) as.numeric(x))
colnames(raw_counts) <- sapply(colnames(raw_counts), function(x) strsplit(x, '/')[[1]][3])
### Sequencing machines add batch and lane info that are not necessary but now don't match covariate file IDs.
seq_ids <- colnames(raw_counts)
colnames(raw_counts) <- alter_seq_ids(colnames(raw_counts))
### If colnames are Sample_Name, change it to Sample_ID
if(all(colnames(raw_counts) %in% covs$Sample_Name)){
  key <- covs$Sample_ID 
  names(key) <- covs$Sample_Name
  colnames(raw_counts) <- key[colnames(raw_counts)]
}

### Account for the annoying case that either of the ID sources has unncessary letter prefaces.
### Returns list of the adjust IDs, if removing some letters from either makes the IDs match exactly when they didn't otherwise.
new_ids <- get_exact_match(colnames(raw_counts), covs$Sample_ID)
covs$ID <- covs$Sample_ID
if(!all(colnames(raw_counts) == new_ids[[1]])){
  colnames(raw_counts) <- new_ids[[1]]
}
if(!all(covs$Sample_ID == new_ids[[2]])){
  covs$ID <- new_ids[[2]]
}
### Match order
row.names(covs) <- covs$ID
covs <- covs[colnames(raw_counts), ]
covs$Seq_ID <- seq_ids

### Add some columns to covs from preprocessing info
covs$N_counts <- colSums(raw_counts)

### Find all STAR log files and read them
star_files <- list.files(path = star_path, pattern = '*/*Log.final.out', recursive = T, full.names = T)
for(star_file in star_files){
  seq_id <- gsub('_Log.final.out', '', tail(strsplit(star_file, '/')[[1]], n = 1))
  star <- read_star_summary(star_file)
  covs[which(covs$Seq_ID == seq_id), 'Aligned_input'] <- star['Number of input reads', 'value']
  covs[which(covs$Seq_ID == seq_id), 'Aligned_uniquely'] <- star['Uniquely mapped reads number', 'value']
  covs[which(covs$Seq_ID == seq_id), 'Aligned_multi'] <- star['Number of reads mapped to multiple loci', 'value']
}

write.table(covs, file = out_file, quote = F, sep = ',', row.names = F, col.names = T)
