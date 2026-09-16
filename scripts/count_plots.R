
library('qqman')
library('stats')
library('ggplot2')
library('methylKit')
library('genomation')
library('GenomicRanges')
library('ggrepel')
library('tools')

if(interactive()){
  project <- '2021612_finch'
  qc_name <- '2024-04-21'
  meth_type <- 'hmC'
  subset_name <- 'none'
  #region <- 'genes'
  region <- 'promoters'
  region_desc <- region
  if(region == 'genes'){
    region_desc <- 'gene_length'
  }
  covs_in <- paste0('/home/cwake/projects/', project, 
                    '/RRBS/SampleSheets/Sample_sheet_2024.csv')
  annot_full <- '/home/cwake/resources/genomes/tguttata/bTaeGut1_v1.p/Annotation/Taeniopygia_guttata.bTaeGut1_v1.p.111.gtf'
  bulk_dir <- '/home/cwake/projects/2021612_finch/results/2024-03-08/'
  
  norm_file <- paste0('/home/cwake/projects/', project, '/RRBS/results/',
                    qc_name, '/methylKit/methylKit_norm_', meth_type, '.RDS')
  adj_file <- paste0('/home/cwake/projects/', project, '/RRBS/results/',
                    qc_name, '/methylKit/norm_batch_', meth_type, '.RDS')
  dm_file <- paste0('/home/cwake/projects/', project, '/RRBS/results/',
                     qc_name, '/methylKit/DM_', meth_type, '_subset-', subset_name, '_region-', region, '.RData')
  
  pdf_file <- paste0('/home/cwake/projects/', project, '/RRBS/results/',
                     qc_name, '/Boxplots_', region_desc, '_', meth_type, '.pdf')
  gtf_rds <- paste0('/home/cwake/projects/', project, '/RRBS/data/gtf.RDS')
  
}else{
  args = commandArgs(trailingOnly=TRUE)
  covs_in <- args[1]
  in_file <- args[2]
  meth_type <- args[3]
  pdf_file <- args[4]
}

remove_dup <- function(proms, dup_loc, step = 1){
  print(dup_loc)
  ### Get reduced proms
  is <- which(paste0(proms@seqnames, '-', proms@ranges@start, '-', proms@ranges@width) == dup_loc)
  is_og <- is
  #d <- proms[which(paste0(proms@seqnames, '-', proms@ranges@start, '-', proms@ranges@width) == dup_loc)]
  len <- length(proms[is]$source)
  biotype_priority <- c('protein_coding', 'lncRNA', 'pseudogene')
  ### If there is more than 1 biotype
  if(length(unique(proms[is]$gene_biotype)) > 1){
    matches <- sapply(proms[is]$gene_biotype, function(x) grep(x, biotype_priority))
    is <- is[which(matches == min(matches))]
    if(length(is) != length(is_og)){
      print('biotype')
    }
  }
  ### Remove those without gene_name
  if(length(proms[is]$source) != 1){
    gene_names <- proms[is]$gene_name
    ### If they aren't both NA, remove those that are
    if(sum(!is.na(gene_names)) != 0){
      is <- is[which(!is.na(proms[is]$gene_name))]
    }
    if(length(is) != length(is_og)){
      print('gene_name')
    }
  }
  if(length(proms[is]$source) != 1){
    if(length(unique(proms[is]$gene_version)) > 1){
      is <- is[which(proms[is]$gene_version == min(proms[is]$gene_version))]
      print('gene_version')
    }
  }
  ### If still there are more than 1, just remove all but the first
  if(length(proms[is]$source) != 1){
    is <- is[1]
    print('No reason at all')
  }
  is_to_remove <- is_og[which(!(is_og %in% is))]
  return(is_to_remove)
}

load(dm_file)
dm <- getData(DM)

print('Reading gtf file')
print(gtf_rds)
if(gtf_rds == '' | is.na(gtf_rds)){
  gtf <- NA
} else{
  if(grepl('\\.gtf', gtf_rds)){
    gtf <- read_gtf(gtf_rds, 'gene', c('gene_name', 'gene_id')) 
  } else{
    gtf <- readRDS(gtf_rds)
  }
  row.names(gtf) <- gtf$gene_id
}
id_to_name <- gtf$gene_name
names(id_to_name) <- gtf$gene_id

### Define bed file options
age <- c('All', 'D6', 'D18', 'D35', 'D90')
#strats <- c(annot_full, paste0(bulk_dir, 'Treatment_Euth_Age-', age, '_DESeq2_sig.bed'))
strats <- c(annot_full, paste0(bulk_dir, 'Treatment_Euth_Age-', age, '_DESeq2_sig.gtf'))
names(strats) <- c('none', age)
print(strats)
print(subset_name)
file.exists(strats)
annot_file <- strats[subset_name]
print(annot_file)

### Read annotation file
### a few tguttata linc RNAs to remove 
linc_reps <- c('ENSTGUG00000024484', 'ENSTGUG00000019743', 'ENSTGUG00000019613', 'ENSTGUG00000018805',
               'ENSTGUG00000026654', 'ENSTGUG00000023892', 'ENSTGUG00000022680', 'ENSTGUG00000025762')
if(file_ext(annot_file) == 'bed'){
  ### GRangesList. exons, introns, promoters, TSSes
  annodat <- readTranscriptFeatures(annot_file)
  #annodat <- readTranscriptFeatures(strats[2])
} else if(file_ext(annot_file) %in% c('gtf', 'gff')){
  ### GRanges object
  genes <- gffToGRanges(annot_file, filter = 'gene')
  to_remove <- which(genes$gene_id %in% linc_reps & genes@strand == '-')
  if(length(to_remove) > 1){
    genes <- genes[-to_remove] 
  }
  ### Default values are upstream = 200, downstream = 200
  proms <- promoters(genes)
  prom_locs <- paste0(proms@seqnames, '-', proms@ranges@start, '-', proms@ranges@width)
  ### If some promters are identical
  if(length(prom_locs) != length(unique(prom_locs))){
    ### Location of those with duplicates
    dup_locs <- names(table(prom_locs))[which(table(prom_locs) > 1)]
    dups <- proms[which(prom_locs %in% dup_locs)]
    ### Return index (in prom) to be removed to make loc unique
    to_remove <- unlist(lapply(dup_locs, function(dup_loc) remove_dup(proms, dup_loc)))
    to_keep <- 1:length(proms$source)
    to_keep <- to_keep[!(to_keep %in% to_remove)]
    proms <- proms[to_keep]
    
    prom_locs <- paste0(proms@seqnames, '-', proms@ranges@start, '-', proms@ranges@width)
    ### Should now be true
    length(prom_locs) == length(unique(prom_locs))
  }

  ### GRangesList
  annodat <- GRangesList(list('genes' = genes,  'promoters' = proms))
}


### Create keys for converting between location and gene_id
an <- as.data.frame(annodat[[region]]@ranges)
### Create key to translate from gene_id to location
id_to_loc <- paste0(annodat[[region]]@seqnames, '-', an$start, '-', an$end)
names(id_to_loc) <- annodat[[region]]$gene_id
### Create key to translate from location to gene_id
loc_to_id <-  annodat[[region]]$gene_id
names(loc_to_id) <- id_to_loc

### Read covariates
covs <- read.csv(covs_in, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
covs <- covs[which(covs$OxBS_or_Mock == 'ox'), ]
row.names(covs) <- covs$Subject_ID

### Assuming format is an RDS file of a methylkit methylbase object
norm <- readRDS(norm_file)
samples <- norm@sample.ids
annot <- covs[samples, c('Treatment', 'Euth_Age', 'Sex', 'Cohort', 'Pool_Plate', 'Library_prep_plate', 'Date_Sequenced', 'FC_ID')]

### Make regional. options are 'site', 'gene', or elements of gtf GRangesList
### ('promoters', 'exons', 'introns', 'TSSes')
if(region %in% c('genes', 'promoters', 'exons', 'introns', 'TSSes')){
  ##### Get counts in 'meth' per region in 'regions'
  ### Save as a database. 
  #meth <- regionCounts(object = meth, regions = annodat[[region]],
  #                     suffix = paste0(names(annot_file), '_', region))
  ### Save as an R object
  meth <- regionCounts(object = norm, regions = annodat[[region]], save.db = F)
  ### Recover gene_id based on their position, as the gene_id is not carried over by regionCounts
  meth@row.names <- loc_to_id[paste0(meth$chr, '-', meth$start, '-', meth$end)]
  meth$gene_names <- id_to_name[meth@row.names]
}

dat <- getData(meth)

coverage <- dat[, colnames(dat)[grepl('coverage', colnames(dat))]]
colnames(coverage) <- samples
numCs <- dat[, colnames(dat)[grepl('numCs', colnames(dat))]]
colnames(numCs) <- samples
dat <- as.data.frame(numCs / coverage)

genes <- c('ENSTGUG00000006498', 'ENSTGUG00000017231', 'ENSTGUG00000000163', 'ENSTGUG00000007140', 'ENSTGUG00000003073')
pdf(pdf_file)
### Boxplots of normalized counts, with DE pvalue and logFC in title
for(id in genes){
  name <-  gtf[id, 'gene_name']
  plot_title <- ifelse(is.na(name), id, name)

  plot_dat <- as.data.frame(covs)
  
  plot_dat[, 'FractionC'] <- as.numeric(dat[id, row.names(plot_dat)])
  plot_dat[, 'Coverage'] <- as.numeric(coverage[id, row.names(plot_dat)])
  plot_dat[, 'Euth_Age'] <- factor(plot_dat[, 'Euth_Age'], levels = c('HATCH', 'FLED', 'IND', 'SM'))
  adjp <- signif(dm[id, 'qvalue'], digits = 3)
  mdiff <- signif(dm[id, 'meth.diff'], digits = 3)
  p1 <- ggplot(plot_dat, aes(x = Euth_Age, y = FractionC)) + 
    labs(title = paste0(plot_title, ' (', meth_type, ', ', gsub('_', ' ', region_desc), ')'),
         subtitle = paste0('MethylKit DM padj: ', adjp, ', meth_diff: ', mdiff)) + xlab('Euth_Age') +
    geom_boxplot(outlier.shape = NA)
  ### Show all points if there are few enough
  if(dim(plot_dat)[1] < 250){
    p1 <- p1 + geom_jitter(color = "red", size = 1, alpha = 0.9)
  }
  print(p1)
  p2 <- ggplot(plot_dat, aes(x = Euth_Age, y = Coverage)) + 
    labs(title = paste0(plot_title, ' (', meth_type, ',', region, ')'),
         subtitle = paste0('MethylKit DM padj: ', adjp, ', meth_diff: ', mdiff)) + xlab('Euth_Age') +
    geom_boxplot(outlier.shape = NA)
  ### Show all points if there are few enough
  if(dim(plot_dat)[1] < 250){
    p2 <- p2 + geom_jitter(color = "red", size = 1, alpha = 0.9)
  }
  print(p2)
}
dev.off()
