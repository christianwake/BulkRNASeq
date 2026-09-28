#library('sys')
library('DESeq2')
library('ggplot2')
library('stringr')
library('GenomicRanges')
library('dplyr')
library('readr')
library('effsize')
library('data.table')
library('biomaRt')
library('fgsea')
library('GSEABase')
library('RColorBrewer')
library('viridis')
library('VennDiagram')
#library('readxl')
library('WriteXLS')
library('circlize')
library('ComplexHeatmap')
library('EnhancedVolcano')

source('/home/cwake/snakemakes/Utility_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')

if(interactive()){
  # project <- '2021612_finch'
  # test <- 'Treatment'
  # strats <- c('Euth_Age-D6', 'Euth_Age-D18','Euth_Age-D35', 'Euth_Age-D90')[1]
  # species <- 'tguttata'
  #custom_set <-  '/home/cwake/projects/2021612_finch/Custom_sets_2022-09-13.xlsx'
  
  project <- '2022612_Petrovas'
  qc_name <- 'QC3'
  test <- 'Strain'
  strats <- c('Cell_type-Tfh', 'Cell_type-pre_Tfh')
  species <- 'mmulatta'
  
  de_files <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/', strats, '/DESeq2_results.txt')
  fgsea_files <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/', strats, '/fgsea_results.txt')
  count_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/finalCounts.txt')
  covs_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_filter.csv')
  gmt_file <- '/home/cwake/resources/gene_sets/c2.cp.v7.2.symbols.gmt'
  out_pdf <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/Heatmaps.pdf')
  pdf_dir <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/genes/')
  custom_set <- ''
}else{
  args = commandArgs(trailingOnly=TRUE)
  res_file <- args[1]
  count_in <- args[2]
  gmt_file <- args[3]
  species <- args[4]
  fgsea_tsv <- args[5]
  out_pdf <- args[6]
  #out_pdf <- args[6]  
}
pthresh <- 0.05
greeks <- c('Α', 'Β', 'Ε', 'Γ', 'Κ')
strat_title <- strsplit(strats, '-')[[1]][1]
### Read covariates
covs <- read.csv(covs_in, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
row.names(covs) <- covs$ID
### Filter samples from covs and counts, and the equivalent in dds
covs <- covs[which(covs$filter == '0'),]
#norm_counts <- norm_counts[, colnames(norm_counts)[which(colnames(norm_counts) %in% row.names(covs))]]

### Read DE results
de_list <- lapply(de_files, function(de_file) read.table(de_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
strats <- gsub(paste0(strat_title, '-'), '', strats)
names(de_list) <- strats
de_genes <- unique(unlist(lapply(de_list, function(x) row.names(x))))
### ID row.names as a column so resilient to rbind
for(i in 1:length(de_list)){
  de_list[[i]]$ID <- row.names(de_list[[i]])
}
de_dat <- as.data.frame(rbindlist(de_list, idcol = strat_title))
de_dat <- de_dat[which(de_dat$padj < 0.05), ]
de_dat$direction <- ifelse(sign(de_dat$log2FoldChange) == 1, '+', '-')

### Confirm direction with AICDA ENSTGUG00000011909
### Read counts
norm_counts <- read.csv(count_in, stringsAsFactors=F, header = T,
                        check.names = F,
                        colClasses = 'character',
                        sep = '\t')
norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))
norm_counts <- norm_counts[, row.names(covs)]

### Manual confirmation of direction
groups <- unique(covs[, test])
d <- c()
for(i in 1:length(row.names(de_dat))){
  g1_mean <- mean(as.numeric(norm_counts[as.character(de_dat[i, 'ID']), covs[which(covs[, test] == groups[1]) , 'Sample_ID']]))
  g2_mean <- mean(as.numeric(norm_counts[as.character(de_dat[i, 'ID']), covs[which(covs[, test]  == groups[2]) , 'Sample_ID']]))
  d <- c(d, ifelse(g2_mean < g1_mean, '-', '+'))
}

de_dat$d = d
all(de_dat$d == de_dat$direction)

de_dat$direction <- ifelse(sign(de_dat$log2FoldChange) == 1, paste0(groups[2], ' is higher (+)'), paste0(groups[1],' is higher (-)'))
table(de_dat[, strat_title], de_dat$direction)

### Read FGSEA results
fgsea <- lapply(fgsea_files, function(fgsea_file) read.table(fgsea_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
names(fgsea) <- strats
fgsea_sig <- lapply(1:length(fgsea), function(i) fgsea[[i]][which(fgsea[[i]]$padj < pthresh), 'pathway'])
names(fgsea_sig) <- strats

names_key <- de_list[[1]]
### Biomart species conversion
if(species != 'hsapiens'){
  ### BiomaRt, for gene names
  id_key <- run_biomaRt(names_key, to_merge = F, species = species)
  ### (For now) Only keep one-to-one orthologs.
  id_key <- id_key[which(id_key$ortholog_type == 'ortholog_one2one'),]
  id_key[id_key == ''] <- NA
  row.names(id_key) <- id_key$ens_short
  colSums(is.na(id_key))
  ### Prioritize Gene names from biomart human orthologs
  names_key[row.names(id_key), 'gene_name'] <- id_key$gene_name
  ### But when that is NA, use gene name from the other species' gtf
  names_key[which(is.na(names_key$gene_name) & !is.na(names_key$gene_name_gtf)), 'gene_name'] <- names_key[which(is.na(names_key$gene_name) & !is.na(names_key$gene_name_gtf)), 'gene_name_gtf']
  names_key[, 'gene_name_plot'] <- names_key$gene_name
  names_key[which(is.na(names_key$gene_name)), 'gene_name_plot'] <- row.names(names_key[which(is.na(names_key$gene_name)), ])
  ### Adjust it if there are duplicates in 'gene_name_plot'
  dups <- names(table(names_key$gene_name_plot))[which(table(names_key$gene_name_plot) > 1)]
  names_key[which(names_key$gene_name_plot %in% dups), 'gene_name_plot'] <- 
    paste0(names_key[which(names_key$gene_name_plot %in% dups), 'gene_name_plot'], '(', row.names(names_key[which(names_key$gene_name_plot %in% dups), ]), ')')
} else{
  names_key$gene_name <- names_key$gene_name_gtf
}
names_key <- names_key[, c('gene_name', 'gene_name_gtf', 'gene_name_plot')]

### Set of genes from the analysis (originally probably ensembl)
goal_set <- names_key$gene_name_gtf[!is.na(names_key$gene_name_gtf)]

if(custom_set != ''){
  ### Read excel file in format - sheet per main category, column per subcatory, rows of variable lengths, genes
  custom <- read_custom_set(custom_set)
  ### For now, in cases of '/', just select the first name even though this indicates aliases already
  custom <- sapply(custom, function(x) strsplit(x, '/')[[1]][1])
  custom <- c(custom, c('p63','p73'))
  ### Returns list of 2 dataframes, output from 'custom' as the 'gene_name' or 'custom' as the 'alias' 
  biomart <- biomaRt_aliases(unique(basic_text_proc(custom, goal_set)))
  ### Fraction of custom set in the goal_set BEFORE matching
  sum(custom %in% goal_set)/length(custom)
  ### Make dataframe 
  mat <- match_sets(custom, goal_set, biomart)
  ### Fraction of custom set in the goal_set AFTER matching
  sum(!is.na(mat$match))/length(custom)
  
  csv_file <- paste0('/home/cwake/projects/', project, '/Custom_sets_matches.csv')
  xls_file <- paste0('/home/cwake/projects/', project, '/Custom_sets_matches.xls')
  write.table(mat, csv_file, quote = F, row.names = F)
  WriteXLS(mat, ExcelFileName = xls_file)
}

##### Choose a set of genes
### From C2CP? It was used in fgsea
DB <- gmtPathways(gmt_file)

pdf(out_pdf)
### Input GMT pathways
geneset_names <- c('REACTOME_MUSCLE_CONTRACTION', 'KEGG_RIBOSOME')
#geneset_names <- c('REACTOME_DNA_METHYLATION', "WP_METHYLATION_PATHWAYS", 'REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA', 'BIOCARTA_BCR_PATHWAY', 'REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE')
for(geneset_name in geneset_names){
  genes <- DB[[geneset_name]]
  genes <- genes[which(genes %in% names_key$gene_name_plot)]
  heatmap_title <- paste0(geneset_name, '\n(', length(genes), ' of ', length(DB[[geneset_name]]), ' genes)')
  names(genes) <- sapply(genes, function(x) row.names(names_key[which(names_key$gene_name_plot == x),]))
  print(length(genes) / length(DB[[geneset_name]]))
  ### fgsea
  if(geneset_name %in% fgsea[[1]]$pathway){
    ps <- sapply(names(de_list), function(i) signif(fgsea[[i]][which(fgsea[[i]]$pathway == geneset_name), 'padj'], digits = 2))
    col_labs <- paste0(names(de_list),'  \n(', ps, ')')  
  } else{
    col_labs <- colnames(de_list)
  }
  ### genes is an array of gene_names (for plotting) with names as IDs as appear in de_list
  my_heatmap(de_list, genes, categories = NA, fontsize = 7, heatmap_title, pheatmap = F, col_clust = F, column_labels = col_labs)
}

dev.off()
### Custom sets
sheets <- excel_sheets(path = custom_set)
### For each sheet in the excel file
for(set_name in sheets){
  set1 <- get_custom_geneset(custom_set, set_name, goal_set, names_key, biomart)
  plot_title <- paste0(set_name, '\n(', length(set1[['genes']]), ' of ',  length(unique(read_custom_set(custom_set, set_name))), ' genes)')
  ### fgsea
  if(set_name %in% fgsea[[1]]$pathway){
    ps <- sapply(names(de_list), function(i) signif(fgsea[[i]][which(fgsea[[i]]$pathway == geneset_name), 'padj'], digits = 2))
    col_labs <- paste0(names(de_list),'  \n(', ps, ')')  
  } else{
    col_labs <- colnames(de_list)
  }
  my_heatmap(de_list, set1[['genes']], set1[['categories']], 7, plot_title, F, col_clust = F, column_labels = col_labs)
}
dev.off()

### Volcano
pthresh <- 0.05
FCthresh <- 2
pdf('/home/cwake/projects/2021612_finch/results/Volcanos.pdf')
for(strat in names(de_list)){
  res <- de_list[[strat]]
  res[row.names(names_key), 'gene_name_plot'] <- names_key$gene_name_plot
  selectLab <- res[which(res$padj <= pthresh & abs(res$log2FoldChange) >= FCthresh & !grepl('^ENSTGUG', res$gene_name_plot)), 'gene_name_plot']
  p <- EnhancedVolcano(res,
                  lab = res[,'gene_name_plot'],
                  selectLab = selectLab,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = strat,
                  caption = '',
                  subtitle = 'Hg/C DE Volcano Plot',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  pCutoff = pthresh,
                  FCcutoff = FCthresh,
                  pointSize = 2.0,
                  labSize = 2.5,
                  colAlpha = 1,
                  legendPosition = 'top',
                  #boxedLabels = T,
                  drawConnectors = TRUE,
                  widthConnectors = 0.5,
                  lengthConnectors = unit(0.01, "npc"),
                  arrowheads = F
  )
  p <- p + theme(legend.position="none")
  print(p)
}
dev.off()
