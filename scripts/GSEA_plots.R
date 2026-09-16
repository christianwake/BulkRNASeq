library('sys')
library('readr')
library('WriteXLS')
library('EnhancedVolcano')
library('data.table')
library('dplyr')
library('VennDiagram')
library('fgsea')
library('circlize')
library('ComplexHeatmap')

source('/home/cwake/snakemakes/Utility_functions.R')
#source('/home/cwake/snakemakes/sc_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')

if(interactive()){
  project <- '2021612_finch'
  qc_name <- '2024-03-08'
  test <- 'Treatment'
  strats_str <- 'All;Euth_Age-D6;Euth_Age-D18;Euth_Age-D35;Euth_Age-D90'
  pthresh <- '0.05'
  
  pdf_out <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/DE_gene_expression.pdf')
  gmt_file <- '/home/cwake/resources/gene_sets/c2.cp.v7.2.symbols.gmt'
  custom_sets <- ''
  count_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/finalCounts.txt')
  covs_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_filter.csv')
  gtf_file <- paste0('/home/cwake/projects/', project, '/data/gtf.RDS')
  de_files <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/', strsplit(strats_str, ';')[[1]], '/DESeq2_results.txt')
  fgsea_files <-  paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/', strsplit(strats_str, ';')[[1]], '/fgsea_results.txt')
}else{
  args = commandArgs(trailingOnly=TRUE)
  pdf_out <- args[1]
  test <- args[2]
  gmt_file <- args[3]
  custom_sets <- args[4]
  strats_str <- args[5]
  count_in <- args[6]
  covs_in <- args[7]
  gtf_file <- args[8]
  res_files <- args[9:length(args)]
  de_files <- res_files[1:(length(res_files)/2)]
  fgsea_files <- res_files[((length(res_files)/2)+1):length(res_files)]
}
pthresh <- 0.05
print(test)
names(test) <- test
if(test == 'Cell_type'){
  names(test) <- 'Infection status'
}
print(strats_str)

strats <- strsplit(strats_str, ';')[[1]]
strat_title <- strsplit(strats[which(strats != 'All')[1]], '-')[[1]][1]
strats <- gsub(paste0(strat_title, '-'), '', strats)
names(strat_title) <- strat_title
if(strat_title == 'Cell_type'){
  names(strat_title) <- 'Infection status'
}

### Read DE results
de_list <- lapply(de_files, function(de_file) read.table(de_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
names(de_list) <- strats
### Read FGSEA results
fgsea <- lapply(fgsea_files, function(fgsea_file) read.table(fgsea_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
names(fgsea) <- strats

### Read covariates
covs <- read.csv(covs_in, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
row.names(covs) <- covs$ID
### Filter samples from covs and counts, and the equivalent in dds
covs <- covs[which(covs$filter == '0'),]
test_desc <- gsub(', ', ' v. ', toString(unique(covs[, test])))

### Confirm direction with AICDA ENSTGUG00000011909
### Read counts
norm_counts <- read.csv(count_in, stringsAsFactors=F, header = T,
                        check.names = F,
                        colClasses = 'character',
                        sep = '\t')
norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))
norm_counts <- norm_counts[, row.names(covs)]

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}
row.names(gtf) <- gtf$gene_id
### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
gtf_cols <- c('gene_name', 'gene_id')
#id_type <- names(which.max(sapply(gtf_cols, function(col) sum(row.names(norm_counts) %in% gtf[, col]))))
id_type <- get_id_name(row.names(norm_counts), gtf)
names_key <- gtf[, c('gene_id', 'gene_name')]
names_key[, 'gene_name_plot'] <- ifelse(is.na(names_key$gene_name), names_key$gene_id, names_key$gene_name)

### Read set of genes
DB <- gmtPathways(gmt_file)

strats <- strats[strats != 'All']
fgsea <- fgsea[strats]
de_list <- de_list[strats]

pdf(pdf_out)

fgsea_sig <- lapply(1:length(fgsea), function(i) fgsea[[i]][which(fgsea[[i]]$padj < pthresh), 'pathway'])
names(fgsea_sig) <- strats

### Checking presence of a particular gene in the sig GSEA results, and their leading edges
gene_of_interest <- 'SLC26A11'
gene_sets <- unlist(fgsea_sig)
gene_sets <- gene_sets[sapply(gene_sets, function(gene_set) gene_of_interest %in% DB[[gene_set]])]
for(gene_set in gene_sets){
  sts <- names(fgsea)[sapply(fgsea, function(fdat) gene_sets[1] %in% fdat$pathway)]
  for(st in sts){
    LE <- strsplit(fgsea[[st]][which(fgsea[[1]]$pathway == gene_set), 'leadingEdgeStr'], ',')[[1]]
    print(length(LE))
    print(gene_of_interest %in% LE)
  }
}

fgsea[[1]]['leadingEdgeStr']

#my_venn(strats, fgsea_sig, paste0('N gene-sets (', test, ' GSEA adjp < ', pthresh, ') compared across ', strat_title))
geneset_names <- unique(unlist(fgsea_sig))
p <- pthresh
while(length(geneset_names) > 50){
  p <- p/10
  print(p)
  fgsea_sig <- lapply(1:length(fgsea), function(i) fgsea[[i]][which(fgsea[[i]]$padj < p), 'pathway'])
  names(fgsea_sig) <- strats
  geneset_names <- unique(unlist(fgsea_sig))
}
if(length(geneset_names) == 0){
  p <- p * 10
  print(p)
  fgsea_sig <- lapply(1:length(fgsea), function(i) fgsea[[i]][which(fgsea[[i]]$padj < p), 'pathway'])
  names(fgsea_sig) <- strats
  geneset_names <- unique(unlist(fgsea_sig))
}
my_venn(strats, fgsea_sig, paste0(test_desc, ' GSEA N gene-sets (adjp < ', p, ')\n compared across ', names(strat_title)))

de_ids <- unique(unlist(lapply(names(de_list), function(x) row.names(de_list[[x]]))))
### Remove genes that weren't in all DE analyses
dat <-as.data.frame(lapply(names(de_list), function(x) de_ids %in% row.names(de_list[[x]])))
row.names(dat) <- de_ids          
de_ids <- de_ids[sapply(row.names(dat), function(rn) rowSums(dat[rn, ]) == length(de_list))]
#geneset_names <- 'REACTOME_RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES'
for(geneset_name in geneset_names){
  print(paste0(geneset_name, ' Heatmap'))
  genes <- DB[[geneset_name]]
  genes <- genes[which(genes %in% gtf$gene_name)]
  column_title <- paste0(test_desc)
  row_title <- paste0(geneset_name, '\n(', length(genes), ' of ', length(DB[[geneset_name]]), ' genes)')
  ### Remove the "[1]" if better handled multiple genes per ID
  names(genes) <- sapply(genes, function(x) row.names(names_key[which(names_key$gene_name_plot == x)[1],]))
  ### Only keep those in the DE results (all)
  genes <- genes[which(names(genes) %in% de_ids)]
  print(length(genes) / length(DB[[geneset_name]]))
  ### fgsea
  if(geneset_name %in% fgsea[[1]]$pathway){
    ps <- sapply(names(de_list), function(i) signif(fgsea[[i]][which(fgsea[[i]]$pathway == geneset_name), 'padj'], digits = 2))
    column_labels <- paste0(names(de_list),' (adj p = ', ps, ')')  
  } else{
    column_labels <- names(de_list)
  }
  fontsize <- heatmap_N_to_size(length(genes))
  ### genes is an array of gene_names (for plotting) with names as IDs as appear in de_list
  my_heatmap(de_list, genes, categories = NA, fontsize = fontsize, column_title, row_title, 
             pheatmap = F, col_clust = F, col_rotation = 7, column_labels = column_labels)
}

dev.off()