#library('sys')
library('DESeq2')
library('ggplot2')
library('stringr')
library('GenomicRanges')
library('dplyr')
library('readr')
library('effsize')
library('data.table')
#library('biomaRt')
library('fgsea')
library('GSEABase')
library('readxl')

source('/home/cwake/snakemakes/Utility_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')
print(sessionInfo())

if(interactive()){
  project <- '2021612_finch'
  qc_name <- '2023-12-06'
  custom_sets <- '/home/cwake/projects/2021612_finch/Custom_sets_2022-09-13.xlsx'
  res_file <- paste0('/home/cwake/projects/', project, 
                     '/results/', qc_name, '/Treatment/All/DESeq2_results.txt')
  
  count_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/finalCounts.txt')
  gtf_file <- paste0('/home/cwake/projects/', project, '/data/gtf.RDS')
  gmt_file <- '/home/cwake/resources/gene_sets/c2.cp.v7.2.symbols.gmt'
  species <- 'mmulatta'
  out_tsv <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/fgsea_Strain_allClusters.tsv')
  out_rds <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/fgsea_Strain_allClusters.RDS')
  out_pdf <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/Strain/Cell_type-pre_Tfh/fgsea_results.pdf')
  #custom_sets <- NA
} else{
  args = commandArgs(trailingOnly=TRUE)
  res_file <- args[1]
  count_in <- args[2]
  gtf_file <- args[3]
  gmt_file <- args[4]
  species <- args[5]
  de_test <- args[6]
  out_tsv <- args[7]
  out_rds <- args[8]
  out_pdf <- args[9]
  custom_sets <- args[10]
}

minSize <- 5

### Read counts
norm_counts <- read.csv(count_in, stringsAsFactors=F, header = T,
                        check.names = F,
                        colClasses = 'character',
                        sep = '\t')
### Convert to numeric
norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))
### Read DE results
res <- read.table(res_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t')

### Read gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

### If gene_name_gtf doesn't exist, add it and (if species != hsapiens) use biomaRt to find human orthologs.
### Then add gene_name_plot.
res <- expand_gene_names(res, species, gtf)

### Remove those without gene name annotation
res <- res[which(!(is.na(res$gene_name))),] 
dup_names <- names(table(res$gene_name))[which(table(res$gene_name) != 1)]
res <- res[which(!res$gene_name %in% dup_names), ]

### Returns a named vector of gene-level stats. Names should be the same as in 'pathways' input to fgsea
ranks  <- get_rank(res, de_test)

### Gene sets
genesets <- getGmt(gmt_file)
DB <- gmtPathways(gmt_file)
### Uses text processing and biomaRt to best match genes in custom_sets and res, then adds those custom sets to DB for use in fgsea
DB <- custom_gene_set_to_biomart_wrapper(DB, custom_sets, res)

### Genes from all  gene sets
#genes <- unique(unlist(lapply(genesets, function(x) geneIds(x))))
genes <- unique(unlist(DB))

### Mitigate sample source bias by making sure that we only have genes with at least one count
a_ids <- row.names(norm_counts)[which(rowSums(norm_counts) > 0)]
a_names <- res[a_ids, 'gene_name']
a_names <- a_names[which(!is.na(a_names))]
paste0('Of ', length(a_ids), ' ids, ', length(a_names), ' are in ', species, '-human biomaRt, and ', sum(a_names %in% genes), ' of those are in human C2CP.')

#DB <- DB[grepl('BIOCARTA', names(DB))]
# redundancy <- unlist(fgseaRes[which(fgseaRes$pathway %in% c('BIOCARTA_P38MAPK_PATHWAY', 'BIOCARTA_TCR_PATHWAY')), 'leadingEdge'])
# le <- unlist(fgseaRes[, 'leadingEdge'])
# a <- sapply(redundancy, function(x) sum(le == x))
# b <- sapply(unique(le), function(x) sum(le == x))
# mean(a)
# mean(b[!names(b) %in% names(a)])

#db <- lapply(names(DB), function(y) DB[[y]][which(DB[[y]] %in% names(ranks))])
#sapply(x, function(y) length(db[[y]]))

print('Beginning fgsea')
fgseaRes <- fgsea(pathways = DB, stats = ranks, minSize = minSize, maxSize = 500)
#fgseaRes <- fgseaRes[, 1:7]
cpdat <- fgseaRes[which(fgseaRes$padj < 0.05),]
cp <- collapsePathways(cpdat, pathways = DB, stats = ranks)
collapsed <- cp[['parentPathways']][which(!is.na(cp[['parentPathways']]))]

fgseaRes <- as.data.frame(fgseaRes)
fgseaRes <- fgseaRes[order(fgseaRes$pval), ]
### Convert the leadingEdge list of gene names into a comma delimited string, in order to print. 
fgseaRes[, 'leadingEdgeStr'] <- apply(fgseaRes, 1, function(x) gsub(', ', ',', toString(x['leadingEdge'][[1]])))

### Add a column with collapsePathways info
row.names(fgseaRes) <- fgseaRes$pathway
fgseaRes[, 'pathway_parent'] <- NA
fgseaRes[cp[['mainPathways']], 'pathway_parent'] <- "Parent"
fgseaRes[names(collapsed), 'pathway_parent'] <- collapsed

### Print significant GSEA results to a file
write.table(fgseaRes[,c("pathway", 'pathway_parent', "pval","padj","ES","NES","size","leadingEdgeStr")], 
            out_tsv, sep = '\t', row.names = F, quote = F)
saveRDS(fgseaRes, file = out_rds)

#fgseaRes <- fgseaRes[, which(colnames(fgseaRes) != 'leadingEdge')]
#WriteXLS(ExcelFileName = out_xls, x = 'fgseaRes', row.names = F, col.names = T)

### Plot of the most significant result
top <- fgseaRes[order(fgseaRes$pval, decreasing = F),][2, 'pathway'][[1]]
ptext <- as.character(signif(fgseaRes[order(fgseaRes$pval, decreasing = F),][2, 'padj'][[1]], digits = 3 ))
pdf(out_pdf)
### Top 10 most significant results
for(i in 1:10){
  ### Plot of the result
  top <- fgseaRes[order(fgseaRes$pval, decreasing = F),][i, 'pathway'][[1]]
  ptext <- as.character(signif(fgseaRes[order(fgseaRes$pval, decreasing = F),][i, 'padj'][[1]], digits = 3 ))
  print(plotEnrichment(DB[[top]], ranks) + labs(title=top, subtitle = paste0(names(res_file),', padj: ',ptext)))
}
dev.off() 

