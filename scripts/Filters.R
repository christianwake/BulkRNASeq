#library('sys')
library('ggplot2')
library('dplyr')
library('genefilter')
library('data.table')
library('scuttle')

source('/home/cwake/snakemakes/Utility_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')

print(Sys.Date())
if(interactive()){
  project <- '2021612_finch'
  qc_name <- '2023-12-06'
  count_file <- paste0('/home/cwake/projects/', project, '/data/counts/featureCounts.txt')
  covs_file <- paste0('/home/cwake/projects/', project, '/data/Covariates_QC_metrics.csv')
  qc_file <- paste0('/home/cwake/projects/', project, '/QC_steps/Sample_and_feature_filters.csv')

  covs_out <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_intermediate.csv')
  out_count <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/data/counts/filteredCounts.txt')
  out_pdf <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/sample_filters.pdf')
  out_txt <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/sample_filters.txt')
  out_txt2 <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/feature_filters.txt')
  chr_pdf <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/chr_dist.pdf')
  
  test <- 'Treatment'
  strat <- 'Euth_Age: D6,D18,D35,D90'
  adjust <- 'Sex,Euth_Age'
  batch <- 'plate, RNA_Extraction_Batch'
  n_reads <- '<=10'
  frac <- '>=0.5'
  keys <- 'RIN_score;Aligned_input;Aligned_uniquely;PCA'
  values <- '<=6;3(SD);3(SD);2(SD)'

  gtf_file <- '/home/cwake/resources/genomes/tguttata/bTaeGut1_v1.p/Annotation/Taeniopygia_guttata.bTaeGut1_v1.p.104.gtf'
} else{
  args = commandArgs(trailingOnly=TRUE)
  count_file <- args[1]
  covs_file <- args[2]
  gtf_file <- args[3]
  qc_file <- args[4]
  covs_out <- args[5]
  out_count <- args[6]
  out_pdf <- args[7]
  out_txt <- args[8]
  out_txt2 <- args[9]
  chr_pdf <- args[10]
  
  test <- args[11]
  strat <- args[12]
  adjust <- args[13]
  batch <- args[14]
}

print('Reading filter file')
filters <- read.table(qc_file, header = T, sep = ',')

### If downample is set to True, any upper limit sample filters for 'N_counts' or  'Aligned_input' will be instead be downsampled to the size of the largest, un-downsampled sample.
### Look for a downsample row and save that information in an object
downsample <- toupper(filters[which(filters$type == 'sample' & tolower(filters$feature) == 'downsample'), 'value'])
downsample <- (downsample == 'TRUE')
### Remove downsample row if it exists
filters <- filters[which(!(filters$type == 'sample' & tolower(filters$feature) == 'downsample')),]

### PCA is handled by the next rule, so remove it now
filters <- filters[which(filters$feature != 'PCA'),]
### Sample filters to a named array 
sample_filters <- filters[which(filters$type  == 'sample'), 'value']
names(sample_filters) <- filters[which(filters$type  == 'sample'), 'feature']
### Save in an array the required covariates (those for a statistical test or sample stratification)
tests <- strsplit(test,',')[[1]]
required <- unique(trimws(c(tests, strsplit(strat,':')[[1]][1], strsplit(adjust,',')[[1]])))

batch <- trimws(strsplit(batch,',')[[1]])
adjust <- trimws(strsplit(adjust, ',')[[1]])
all_covs <- c(tests, adjust, batch)
all_covs <- all_covs[which(all_covs != '')]

### Set defaults
### If no version of depth filter is input, add N_reads default
if(!any(c('Aligned_input', 'Aligned_uniquely', 'N_reads', 'N_counts') %in% names(sample_filters))){
  x <- '3(SD),3(SD)'
  names(x) <- 'N_reads'
  sample_filters <- c(sample_filters, x)
}

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

### Read in covariates
covs <- read.csv(covs_file, stringsAsFactors=F, header = T,
                       check.names=F, 
                       colClasses = 'character',
                       sep=',')
row.names(covs) <- covs$ID
### Order, sequencing postscripts and differing subscripts were already accounted for in Compile_QC_metrics.R
colnames(raw_counts) <- covs$ID
covs[, all_covs] <- sapply(all_covs, function(col) gsub(' ', '_', covs[, col]))

### Filter those with missing necessary information
covs <- covs[which(rowSums(is.na(covs[, required])) == 0),]

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
  #id_type <- names(which.max(sapply(gtf_cols, function(col) sum(row.names(raw_counts) %in% gtf[, col]))))
  id_type <- get_id_name(row.names(raw_counts), gtf)
}

### Remove specified samples
if('Sample_Name' %in% names(sample_filters)){
  id <- covs[which(covs$Sample_Name == sample_filters['Sample_Name']), 'ID']
  ### Remove from covs
  covs <- covs[which(covs$Sample_Name != sample_filters['Sample_Name']),]
  ### Remove from raw_counts
  raw_counts <- raw_counts[, which(colnames(raw_counts) != id)]
  print('Removing sample(s) based on input Sample_Name')
}
if('Sample_ID' %in% names(sample_filters)){
  ### Remove from covs
  covs <- covs[which(covs$Sample_Name != sample_filters['Sample_ID']),]
  ### Remove from raw_counts
  raw_counts <- raw_counts[, which(colnames(raw_counts) != sample_filters['Sample_ID'])]
  print('Removing sample(s) based on input Sample_ID')
}

### RNA Integrity Number
rins <- c('RIN', 'RIN_SCORE', 'RNA_INTEGRITY')
rin <- names(sample_filters)[toupper(names(sample_filters)) %in% rins]
if(length(rin) == 0){
  warning('Is no RIN filter entered?')
} else{
  print('Removing samples based on RIN score')
  ### How to handle missing RIN info?
  covs[, rin]<- suppressWarnings(as.numeric(covs[, rin]))
  ### If the fraction of missing RIN values is less than 0.35, replace those with the mean of the present RINs.
  if(sum(is.na(covs[, rin]))/length(covs[, rin]) < 0.35){
    covs[which(is.na(covs[, rin])), rin] <- mean(covs[, rin], na.rm = T)
  }
}

### If 'Sex' is among the required covariates, filter those that are not clearly 'Male' or 'Female', e.g. NA or UNK
filter_sexless <- T
sex <- required[which(toupper(required) %in% c("SEX","GENDER"))]
print('Checking sex')
if(length(sex) > 0){
  sexes <- names(table(covs[, sex]))
  if(length(sexes) > 2){
    odd <- sexes[which(!toupper(sexes) %in% c('M','F','MALE','FEMALE','0','1','2'))]
    ### If the remaining are binary
    if(length(sexes[which(!sexes %in% odd)]) == 2){
      if(filter_sexless){
        covs <- covs[which(!covs[, sex] %in% odd),]
        raw_counts <- raw_counts[, covs$ID]
      }
    }
  }
}

### Make numeric
cols <- c('RIN_score', 'N_counts', 'Aligned_input', 'Aligned_uniquely', 'Aligned_multi')
cols <- cols[which(cols %in% colnames(covs))]
covs <- covs[, colnames(covs)]
for(c in cols){
  covs[, c] <- as.numeric(covs[, c])
}
covs[, 'nFeature'] <- colSums(raw_counts > 0) 

print('mito')
### Determine which genes are mitochondrial by gtf if it is input. Otherwise, by gene name pattern '^MT'.
if(any(c('MT', 'mito', 'MT_sum', 'MT_frac', 'MT_prop') %in% filters$feature)){
  if(!is.na(gtf)){
    covs[,'MT_sum'] <- sum_chromosome(raw_counts, chr ='MT', gtf, id_type)
  } else {
    mts <- row.names(raw_counts)[grepl('^MT-', row.names(raw_counts))]
    covs['MT_sum'] <- colSums(raw_counts[mts, ])
  }
}

# feats <- rowSums(raw_counts)
# feats <- feats[order(feats, decreasing = T)]
# gtf[names(feats[1:100]), 'gene_name']

#plot(raw_counts[gtf[which(gtf$gene_name == 'RPL13'), 'gene_id'], ], as.numeric(covs[,'N_counts']))
print('Chromosome distribution')
chrs <- unique(gtf$chr)
chrmax <- suppressWarnings(max(as.numeric(chrs), na.rm = T))
chr_dat <- lapply(chrs, function(chr) as.data.frame(t(as.data.frame(sum_chromosome(raw_counts, chr = chr, gtf, id_type)))))
chr_dat <- as.data.frame(rbindlist(chr_dat))
row.names(chr_dat) <- chrs
chr_dat <- chr_dat[c(as.character(1:chrmax), 'X', 'Y', 'MT'),]
pdf(chr_pdf)
for(col in colnames(chr_dat)){
  df <- chr_dat[, col, drop = F]
  colnames(df) <- c('counts')
  df$chr <- row.names(df)
  p <- ggplot(data= df , aes(x = chr, y = counts)) + ggtitle(col) +
  geom_bar(stat = "identity") + scale_x_discrete(limits=row.names(df))
  print(p)
}
dev.off()

#####
##### Feature filters
#####
print('Beginning feature filters')
n_reads <- filters[which(filters$type == 'feature' & filters$feature == 'counts'), 'value']
frac <- filters[which(filters$type == 'feature' & filters$feature == 'fraction_samples'), 'value']

if(n_reads == '' & frac == ''){
  frac <- '>=0.9'
}

exclude_gene_names <- c()
feat_txt <- filters[which(filters$type == 'feature'),]
feat_txt$N <- NA

print('N_reads')
small_features <- c()
if(n_reads != ''){
  ### Input can be 1 or 2 (comma-delimited) values. If 2, first is assumed to be lower-bound, second to be upper-bound.
  ### Values can be numerical and be followed by a '(SD)', can be preceded by a <, >, =< or => (if only 1 input value).
  thresh <- threshold_string_dat(covs, n_reads, 'N_reads')
  if(grepl('=', n_reads)){
    small_features <- rownames(raw_counts)[rowSums(raw_counts) <= thresh[1]]
  } else{
    small_features <- rownames(raw_counts)[rowSums(raw_counts) < thresh[1]]
  }
  ### Boxplot the distribution
  dat <- as.data.frame(rowSums(raw_counts))
  colnames(dat) <- 'N_reads'
  ggplot(dat, aes(x = log(N_reads))) + 
    geom_histogram(fill="slateblue")+ 
    geom_vline(xintercept = log(thresh[1]), col = 'red') +
    ggtitle(paste0('Distr. of N reads across features\n Threshold.', n_reads))
  
  exclude_gene_names <- c(exclude_gene_names, small_features)
  feat_txt[which(feat_txt$feature == 'counts'), 'N'] <- length(small_features) 
}

print('fraction 0')
sparse_features <- c()
if(frac != ''){
  ### frac is fraction of samples with 0 counts for that transcript. e.g '>=0.2' will select for filtering transcripts with 0 counts in more than 20% of samples
  sparse_features <- Gene_filter(raw_counts, as.numeric(gsub('>', '', gsub('<', '', gsub('=', '', frac)))), or_equals = grepl('=', frac), return_features = T)
  ### Boxplot the distribution
  dat <- as.data.frame(rowSums(raw_counts == 0))
  colnames(dat) <- 'N_empty'
  ggplot(dat, aes(x = N_empty)) + 
    geom_histogram(fill="slateblue")+ 
    geom_vline(xintercept = as.numeric(gsub('<', '', gsub('=', '', frac))) * length(colnames(raw_counts)), col = 'red') +
    ggtitle(paste0('Distr. of 0-count samples across features. \n Threshold ', frac, ''))
  
  feat_txt[which(feat_txt$feature == 'fraction_samples'), 'N'] <- length(sparse_features) 
  exclude_gene_names <- c(exclude_gene_names, sparse_features)
}

print('variance')
var_thresh <- filters[which(filters$type == 'feature' & filters$feature == 'variance'), 'value']
### Based on variance
if(var_thresh != ''){
  features <- row.names(raw_counts)[which(!row.names(raw_counts) %in% exclude_gene_names)]
  vars <- sapply(features, function(rn) var(as.numeric(raw_counts[rn, ])))
  #vars <- vars[order(vars)]
  thresh <- threshold_string_feature(vars, var_thresh)
  if(grepl('=', var_thresh)){
    stable_features <- rownames(vars)[vars <= thresh[1]]
  } else{
    stable_features <- rownames(vars)[vars < thresh[1]]
  }
  feat_txt[which(feat_txt$feature == 'variance'), 'N'] <- length(stable_features) 
  exclude_gene_names <- c(exclude_gene_names, stable_features)
}

print('key word')
### Based on input key words
key_words <- filters[which(filters$type == 'feature' & filters$feature == 'pattern'), 'value']
if(key_words != ''){
  key_words <- strsplit(key_words, ',')[[1]]
  if(id_type == 'gene_id'){
    all_names <-  gtf$gene_name
    genes <- row.names(raw_counts)[sapply(row.names(raw_counts), function(id) any(sapply(key_words, function(y) grepl(pattern = y, x = gtf[id, 'gene_name'], perl = T))))]
  } else{
    all_names <- row.names(raw_counts)
    genes <- all_names[sapply(all_names, function(x) any(sapply(key_words, function(y) grepl(pattern = y, x = x, perl = T))))]
  }
  
  exclude_gene_names <- c(exclude_gene_names, genes)
  feat_txt[which(feat_txt$feature == 'pattern'), 'N'] <- length(genes) 
}

exclude_gene_names <- sort(unique(exclude_gene_names))
feat_txt <- rbind(feat_txt, c('feature', 'total filtered', '', length(unique(exclude_gene_names))))
feat_txt <- rbind(feat_txt, c('feature', 'total remaining', '', length(row.names(raw_counts)) - length(unique(exclude_gene_names))))


write.table(feat_txt, out_txt2, quote = F, sep = '\t', row.names = T, col.names = T)

### Filter counts
raw_counts <- raw_counts[which(!row.names(raw_counts) %in% unique(exclude_gene_names)), ]
covs$N_counts <- colSums(raw_counts)

### Separate flat from SD filter. Couldn't get this right, so they must be entered in the config separately
# spl <- strsplit(sample_filters[1], ',')[[1]]
# for(i in 1:length(sample_filters)){
#   ### If there is more than one entry (upper and lower)
#   if(length(strsplit(sample_filters[i], ',')[[1]]) > 1){
#     ### If upper and lower differ by flat/SD
#     if(length(unique(sapply(strsplit(sample_filters[i], ',')[[1]], function(x) grepl('\\(SD\\)', x)))) != 1){
#       sample_filters <- c(sample_filters, paste0('<=', strsplit(sample_filters[i], ',')[[1]][1]), paste0('>=', strsplit(sample_filters[i], ',')[[1]][2]))
#     }
#   }
# }
# dat <- as.data.frame(t(as.data.frame(lapply(names(sample_filters), function(x) strsplit(sample_filters[x], ',')[[1]]))))
# colnames(dat) <- c('upper', 'lower')
# dat[, 'feature'] <- names(sample_filters)

print('Sample filters')

flat <- sample_filters[!sapply(sample_filters, function(x) grepl('\\(SD\\)', x))]
sds <- sample_filters[sapply(sample_filters, function(x) grepl('\\(SD\\)', x))]


pdf(out_pdf)
plot(covs$nFeature, covs$Aligned_input, xlab = 'Number of features > 1', ylab= 'N aligned reads')
plot(covs$nFeature, covs$N_counts, xlab = 'Number of features > 1', ylab= 'N counts')

# ggplot(covs, aes(x=nFeature, y=N_counts)) +
#   geom_point() + # Show dots
#   geom_text(label = gsub('S', '', row.names(covs)), nudge_y = 2000000, check_overlap = F) + 
#   geom_hline(yintercept = 1000000, linetype="solid", color = "red", size=0.5)

covs_full <- covs
to_filter <- c()
to_downsample <- c()

txt_flat <- as.data.frame(matrix(ncol = 4, nrow = length(flat)))
colnames(txt_flat) <- c('filter', 'n', 'lower_failed', 'upper_failed')
txt_flat[, 'filter'] <- names(flat)
row.names(txt_flat) <- names(flat)
### Apply flat filters
print('Flat filters first')
if(length(flat) > 0){
  for(i in 1:length(flat)){
    thresh <- threshold_string_dat(covs, flat[i], names(flat)[i])
    if(grepl('=', flat[i])){
      lower <- covs[which(covs[, names(flat)[i]] <= thresh[1]), 'ID']
      upper <- covs[which(covs[, names(flat)[i]] >= thresh[2]), 'ID']
    } else{
      lower <- covs[which(covs[, names(flat)[i]] < thresh[1]), 'ID']
      upper <- covs[which(covs[, names(flat)[i]] > thresh[2]), 'ID']
    }
    ### Add N filtered to txt object
    txt_flat[names(flat)[i], 'n'] <- length(c(lower, upper))
    txt_flat[names(flat)[i], 'lower_failed'] <- gsub(',', ';', toString(lower))
    txt_flat[names(flat)[i], 'upper_failed'] <- gsub(',', ';', toString(upper))
    ### Boxplot the distribution
    p <- ggplot(covs, aes(y=eval(parse(text=names(flat))))) + 
      geom_boxplot(fill="slateblue", alpha=0.2) + 
      ylab(names(flat))
    ### Add filter threshold lines
    for(t in thresh[!is.na(thresh)]){
      p <- p + geom_hline(yintercept = t, col = 'red')
    }
    print(p)
  }  
}

### Add samples to to_filter (except in case of downsampling instead of filtering)
for(rn in row.names(txt_flat)){
  to_filter <- c(to_filter, strsplit(txt_flat[rn, 'lower_failed'], '; ')[[1]])
  ### Add upper limit samples to be filtered (UNLESS downsampled is set to TRUE and this filter is a size filter)
  if(downsample != T | !(txt_flat[rn, 'filter'] %in% c('N_counts', 'Aligned_input'))){
    to_filter <- c(to_filter, strsplit(txt_flat[rn, 'upper_failed'], '; ')[[1]])
  } else{
    to_downsample <- c(to_downsample, strsplit(txt_flat[rn, 'upper_failed'], '; ')[[1]])
  }
}

### Filter from covs dataframe
covs <- covs[which(!row.names(covs) %in% to_filter), ]

print('Standard deviation filters second')
txt_sd <- as.data.frame(matrix(ncol = 4, nrow = length(sds)))
colnames(txt_sd) <- c('filter', 'n', 'lower_failed', 'upper_failed')
txt_sd[, 'filter'] <- names(sds)
row.names(txt_sd) <- names(sds)

### Apply SD filters
to_filter <- c()
if(length(sds) > 0 ){
  threshes <- as.list(rep(NA, length(sds)))
  names(threshes) <- names(sds)
  for(i in 1:length(sds)){
    thresh <- threshold_string_dat(covs, sds[i], names(sds)[i])
    threshes[[names(sds)[i]]] <- thresh
    if(grepl('=', sds[i])){
      lower <- covs[which(covs[, names(sds)[i]] <= thresh[1]), 'ID']
      upper <- covs[which(covs[, names(sds)[i]] >= thresh[2]), 'ID']
    } else{
      lower <- covs[which(covs[, names(sds)[i]] < thresh[1]), 'ID']
      upper <- covs[which(covs[, names(sds)[i]] > thresh[2]), 'ID']
    }
    ### Add N filtered to txt object
    txt_sd[names(sds)[i], 'n'] <- length(c(lower, upper))
    txt_sd[names(sds)[i], 'lower_failed'] <- gsub(',', ';', toString(lower))
    txt_sd[names(sds)[i], 'upper_failed'] <- gsub(',', ';', toString(upper))
    ### Boxplot the distribution
    p <- ggplot(covs, aes(y=eval(parse(text=names(sds)[i])))) + 
      geom_boxplot(fill="slateblue", alpha=0.2) + 
      ylab(names(sds)[i])
    ### Add filter threshold lines
    for(t in thresh[!is.na(thresh)]){
      p <- p + geom_hline(yintercept = t, col = 'red')
    }
    print(p)
  }
}

### Add samples to to_filter (except in case of downsampling instead of filtering)
for(rn in row.names(txt_sd)){
  to_filter <- c(to_filter, strsplit(txt_sd[rn, 'lower_failed'], '; ')[[1]])
  ### Add upper limit samples to be filtered (UNLESS downsampled is sest to TRUE and this filter is a size filter)
  if(downsample != T | !(txt_sd[rn, 'filter'] %in% c('N_counts', 'Aligned_input'))){
    to_filter <- unique(c(to_filter, strsplit(txt_sd[rn, 'upper_failed'], '; ')[[1]]))
  } else{
    to_downsample <- unique(c(to_downsample, strsplit(txt_sd[rn, 'upper_failed'], '; ')[[1]]))
  }
}

txt <- rbind(txt_flat, txt_sd)
### Filter from covs dataframe
covs <- covs[which(!row.names(covs) %in% to_filter), ]
txt['Remaining', ] <- c('Remaining', dim(covs)[1], NA, NA)
### Remove those from raw_counts
raw_counts <- raw_counts[, which(colnames(raw_counts) %in% row.names(covs))]

print('Down sampling read count upper limit samples (if specified)')
### Apply downsampling (by N raw reads)
if(length(to_downsample) > 0){
  ### Count sums of eachsample
  csums <- colSums(raw_counts)
  ### Order it by read count
  csums <- csums[order(csums, decreasing = T)]
  to_downsample <- names(csums)[which(names(csums) %in% to_downsample)]
  ### Maximum read count that won't be downsampled (for reference)
  ref1 <- names(csums[!(names(csums) %in% to_downsample)])[1]
  
  ### Get the proportion to be used for downsampling to bring each sample down to the read counts of ref1
  ds_props <- sapply(to_downsample, function(s) 1 - (csums[s] - csums[ref1])/csums[s])
  names(ds_props) <- to_downsample
  ### Adjust slightly to keep order
  ds_props <- ds_props + seq(from = (length(ds_props) * 0.01), to = 0.01, by = -0.01)
  ### Full array of proportions, for the downsampleMatrix function
  props <- rep(1, length(csums))
  names(props) <- names(csums)
  ### Add downsampling proportion to the array of 1s
  props[to_downsample] <- ds_props
  ### Match order to raw_counts
  props <- props[colnames(raw_counts)]
  ### Downsample
  #raw_counts2 <- downsampleMatrix(raw_counts, prop = props, bycol = T)
  raw_counts <- downsampleMatrix(raw_counts, prop = props, bycol = T)
  raw_counts <- as.matrix(raw_counts)
  new_csums <- colSums(raw_counts)
  new_csums <- new_csums[names(csums)]
  
  c1 <- as.data.frame(csums)
  c1$type <- 'Raw'
  colnames(c1) <- c('counts', 'type')
  c1$sample <- row.names(c1)
  
  c2 <- as.data.frame(new_csums)
  c2$type <- 'Downsampled'
  colnames(c2) <- c('counts', 'type')
  c2$sample <- row.names(c2)
  
  plot_dat <- rbind(c1, c2)
  

  p <- ggplot(plot_dat, aes(x = counts, y = sample)) + 
    geom_line() + theme(axis.text=element_text(size = 6)) + ggtitle('Downsampling') + 
    geom_point(size = 2, aes(color = type)) +
    geom_vline(xintercept = threshes[['N_counts']][2], linetype = 'dashed') +
    geom_vline(xintercept = csums[ref1], linetype = 'solid')
  print(p)
}

### Add a filter column to covs and print, or remove those samples and print
covs_full[, 'filter'] <- ifelse(covs_full$ID %in% covs$ID, 0, 1)
covs_full[, 'downsampled'] <- sapply(covs_full$ID, function(x) ifelse(x %in% to_downsample, ds_props[x], 0))
### Write new covs and summary txt
write.table(txt, out_txt, quote = F, sep = ',', row.names = F, col.names = T)
write.table(covs_full, covs_out, quote = F, sep = ',', row.names = F, col.names = T)
write.table(raw_counts, out_count, quote = F, sep = '\t', row.names = T, col.names = T)

plot(covs$nFeature, covs$Aligned_input, xlab = 'Number of features > 1', ylab= 'N aligned reads')

dev.off()
print(table(covs[,test], covs[, strsplit(strat, ':')[[1]][1]]))
