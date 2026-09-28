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
  out_tsv <- args[6]
  out_rds <- args[7]
  out_pdf <- args[8]
  custom_sets <- args[9]
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
### Run biomaRt to get gene name information from ensembl IDs.
res$ensembl <- row.names(res)
#colnames(res) <- gsub('gene_name', 'gene_name_gtf', colnames(res))

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
  id_type <- get_id_name(row.names(norm_counts), gtf)
}
### Determine whether each is a gene ID or name
all_names <- sapply(row.names(res), function(x) get_id_name(x, gtf))
names(all_names) <- row.names(res)
table(all_names)
### For each gene name or ID, if it is a gene name, use it. If not, get the IDs gene name match from gtf, or just the ID if there is none
res$gene_name <- sapply(1:length(all_names), function(i) ifelse(all_names[i] == 'gene_name', 
                                               names(all_names)[i], 
                                               ifelse(is.na(gtf[names(all_names)[i], 'gene_name']),
                                                      names(all_names)[i],gtf[names(all_names)[i], 'gene_name'])))
res$gene_name_gtf <- sapply(1:length(all_names), function(i) ifelse(all_names[i] == 'gene_name', 
                                                                names(all_names)[i], 
                                                                gtf[names(all_names)[i], 'gene_name']))
id_key <- res[, 'gene_name', drop= F]

# if(species != 'hsapiens'){
#   ### BiomaRt, for gene names
#   id_key <- run_biomaRt(res, to_merge = F, species = species)
#   ### (For now) Only keep one-to-one orthologs.
#   id_key <- id_key[which(id_key$ortholog_type == 'ortholog_one2one'),]
#   id_key[id_key == ''] <- NA
#   row.names(id_key) <- id_key$ens_short
#   colSums(is.na(id_key))
#   ### Prioritize Gene names from biomart human orthologs
#   res[row.names(id_key), 'gene_name'] <- id_key$gene_name
#   ### But when that is NA, use gene name from the other species' gtf
#   res[which(is.na(res$gene_name) & !is.na(res$gene_name_gtf)), 'gene_name'] <- res[which(is.na(res$gene_name) & !is.na(res$gene_name_gtf)), 'gene_name_gtf']
#   res[, 'gene_name_plot'] <- res$gene_name
#   res[which(is.na(res$gene_name)), 'gene_name_plot'] <- row.names(res[which(is.na(res$gene_name)), ])
#   ### Adjust it if there are duplicates in 'gene_name_plot'
#   dups <- names(table(res$gene_name_plot))[which(table(res$gene_name_plot) > 1)]
#   res[which(res$gene_name_plot %in% dups), 'gene_name_plot'] <- 
#     paste0(res[which(res$gene_name_plot %in% dups), 'gene_name_plot'], '(', row.names(res[which(res$gene_name_plot %in% dups), ]), ')')
#   
# } else{
#   res$gene_name <- res$gene_name_gtf
# }

### Remove those without gene name annotation
res <- res[which(!(is.na(res$gene_name))),] 
dup_names <- names(table(res$gene_name))[which(table(res$gene_name) != 1)]
res <- res[which(!res$gene_name %in% dup_names),]
#rnk <- deseqRes2Rnk(res, by = '-log(p)', id = 'gene_name')
#rnk <- rnk[order(rnk$NegLogP, decreasing = T),]
#ranks <- setNames(rnk[, 'NegLogP'], rnk$ID)
rnk <- deseqRes2Rnk(res, by = 'stat', id = 'gene_name')
rnk <- rnk[order(rnk$stat, decreasing = T),]
ranks <- setNames(rnk[, 'stat'], rnk$ID)

### Gene sets
genesets <- getGmt(gmt_file)
DB <- gmtPathways(gmt_file)
### Custom sets
goal_set <- res$gene_name_gtf[!is.na(res$gene_name_gtf)]
greeks <- c('Α', 'Β', 'Ε', 'Γ', 'Κ')

if(custom_sets != '' & !is.na(custom_sets)){
  sheets <- excel_sheets(path = custom_sets)
  ### Do Biomart all at once and save it, because it often has connection issues and will fail, so this simplifies that.
  biomart_file <- gsub('.xlsx', '.RDS', custom_sets)
  if(!file.exists(biomart_file)){
    biomart_custom_geneset(custom_sets, sheets, goal_set, biomart_file)
  }
  ### When biomart is NA, it runs biomaRt_aliases(cust) to get it. 
  ### When biomart is an RDS file, it simply reads it
  cus <- lapply(sheets, function(set_name) get_custom_geneset(
    custom_sets, set_name, goal_set, res, biomart = biomart_file)[['genes']])
  names(cus) <- sheets
  ### Add custom sets to the DB object
  DB <- c(DB, cus)
} 

### Genes from all  gene sets
#genes <- unique(unlist(lapply(genesets, function(x) geneIds(x))))
genes <- unique(unlist(DB))

### Mitigate sample source bias by making sure that we only have genes with at least one count
a_ids <- row.names(norm_counts)[which(rowSums(norm_counts) > 0)]
a_names <- id_key[a_ids, 'gene_name']
a_names <- a_names[which(!is.na(a_names))]

#paste0('Of ', length(a_ids), ' ids, ', length(a_names), ' are in ', species, '-human biomaRt, and ', sum(a_names %in% genes), ' of those are in human C2CP.')

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

