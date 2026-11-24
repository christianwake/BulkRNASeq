#library('sys')
#library('readr')
library('WriteXLS')
library('EnhancedVolcano')
library('data.table')
library('dplyr')
library('VennDiagram')
library('circlize')
library('ComplexHeatmap')
#library('pheatmap')
library('RColorBrewer')

source('/home/cwake/snakemakes/Utility_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')

if(interactive()){
  project <- '2021612_finch'
  qc_name <- '2024-03-08'
  test <- 'Treatment'
  pthresh <- '0.05'
  strats_str <- 'All;Euth_Age-D6;Euth_Age-D18;Euth_Age-D35;Euth_Age-D90'
  
  pdf_out <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, 'DE2.pdf')
  count_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/counts/finalCounts.txt')
  covs_in <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/QC/Covariates_QC_metrics_filter.csv')
  gtf_file <- paste0('/home/cwake/projects/', project, '/data/gtf.RDS')
  de_files <- paste0('/home/cwake/projects/', project, '/results/', qc_name, '/', test, '/', strsplit(strats_str, ';')[[1]], '/DESeq2_results.txt')
  
}else{
  args = commandArgs(trailingOnly=TRUE)
  pdf_out <- args[1]
  test <- args[2]
  strats_str <- args[3]
  count_in <- args[4]
  covs_in <- args[5]
  gtf_file <- args[6]
  de_files <- args[7:length(args)]
}
pthresh <- '0.001'

print(test)
names(test) <- test
if(test == 'Cell_type'){
  names(test) <- 'Infection status'
}

print(strats_str)
print(de_files)
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
de_genes <- unique(unlist(lapply(de_list, function(x) row.names(x))))
### ID row.names as a column so resilient to rbind
for(i in 1:length(de_list)){
  de_list[[i]]$ID <- row.names(de_list[[i]])
}
de_dat <- as.data.frame(rbindlist(de_list, idcol = names(strat_title)))
de_dat <- de_dat[which(de_dat$padj < 0.05), ]
de_dat$direction <- ifelse(sign(de_dat$log2FoldChange) == 1, '+', '-')

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

# strats <- strats[strats != 'All']
# de_list <- de_list[strats]

pdf(pdf_out)
ven <- lapply(de_list, function(x) row.names(x[which(x$padj < as.numeric(pthresh)),]))
my_venn(strats, ven, paste0(test_desc, ' N DEG (adjp < ', pthresh, ')\n compared across ', names(strat_title)))

### Heatmap of LogFC values for multiple DE comparisons (columns)
run_heatmap <- function(de_list, names_key, test_desc, pthresh){
  ### sig IDs in any DE analysis
  de_ids <- unique(unlist(lapply(names(de_list), function(x) row.names(de_list[[x]][which(de_list[[x]][, 'padj'] < pthresh),]))))
  ### Remove genes that weren't in all DE analyses
  dat <- as.data.frame(lapply(names(de_list), function(x) de_ids %in% row.names(de_list[[x]])))
  row.names(dat) <- de_ids          
  de_ids <- de_ids[sapply(row.names(dat), function(rn) rowSums(dat[rn, ]) == length(de_list))]
  genes <- names_key[de_ids, 'gene_name_plot']
  names(genes) <- de_ids
  column_title <- paste0(test_desc)
  row_title <- paste0('DEG adj p < ', pthresh)
  ### Remove the "[1]" if better handled multiple genes per ID
  names(genes) <- sapply(genes, function(x) row.names(names_key[which(names_key$gene_name_plot == x)[1],]))
  ### Column labels
  #column_labels <- paste0('in ', names(de_list))
  ### genes is an array of gene_names (for plotting) with names as IDs as appear in de_list
  fontsize <- heatmap_N_to_size(length(genes))
  p0 <- my_heatmap(de_list, genes, categories = NA, fontsize = fontsize, column_title, row_title, pheatmap = F, 
             col_clust = F)
  return(genes)
}

genes <- run_heatmap(de_list[strats != 'All'], names_key, test_desc, as.numeric(pthresh))
p <- as.numeric(pthresh)
while(length(genes) > 50){
  p <- p/10
  genes <- run_heatmap(de_list[strats != 'All'], names_key, test_desc, p)
}

### Heatmap of Normalized Expression values (genes x subject with subject annotations)
col_clust <- F
col_dend <- F

#plot_dat <- as.matrix(log(norm_counts[names(genes), ] + 1))
plot_dat <- as.matrix(norm_counts[names(genes), ])
row.names(plot_dat) <- genes[row.names(plot_dat)]
### Reorder columns
covs <- covs[order(covs[, strat_title], covs[, test]), ]
plot_dat <- plot_dat[, row.names(covs)]

column_labels <- covs$Sample_Name
### Annotation Colors
cats <- c(unique(covs[, test]), unique(covs[, strat_title]))
col_color <- brewer.pal(n = length(cats), 'Set2')
names(col_color) <- cats
col_color <- list(Function = col_color)
ComplexHeatmap_anno <- columnAnnotation(df = covs[, c(test, strat_title), drop = F], col = col_color)

### genes is an array of gene_names (for plotting) with names as IDs as appear in de_list
fontsize <- heatmap_N_to_size(length(genes))

### Scale
scaled_dat <- t(scale(t(plot_dat)))
Min <- min(plot_dat)
Max <- max(plot_dat)
### Cell colors
ComplexHeatmap_cols <- colorRamp2(c(Min, 0, Max), c("blue", "white", "red"))

p1 <- Heatmap(plot_dat,
             column_title = paste0(test_desc, ' DE genes'),
             row_names_gp = gpar(fontsize = fontsize), ### Font size of row labels
             col = ComplexHeatmap_cols, ### Function to create colors per value,
             name = 'norm. \nexpr.', ### Color legend title
             cluster_columns = col_dend,
             column_dend_reorder = F, ### Do clustering and keep dendrogram, but don't order by the cluster (pheatmap can't do this)
             bottom_annotation = ComplexHeatmap_anno,
             column_labels = column_labels,
             column_names_rot = 60
             
)
draw(p1, show_heatmap_legend = T) 

Min <- min(scaled_dat)
Max <- max(scaled_dat)
### Cell colors
ComplexHeatmap_cols <- colorRamp2(c(Min, 0, Max), c("blue", "white", "red"))

p2 <- Heatmap(scaled_dat,
             column_title = paste0(test_desc, ' DE genes'),
             row_names_gp = gpar(fontsize = fontsize), ### Font size of row labels
             col = ComplexHeatmap_cols, ### Function to create colors per value,
             name = 'scaled\nnorm. \nexpr.', ### Color legend title
             cluster_columns = col_dend,
             column_dend_reorder = F, ### Do clustering and keep dendrogram, but don't order by the cluster (pheatmap can't do this)
             bottom_annotation = ComplexHeatmap_anno,
             column_labels = column_labels,
             column_names_rot = 60
             
)
draw(p2, show_heatmap_legend = T)


### Volcano plots
FCthresh <- 2
pthresh <- as.numeric(pthresh)
for(strat in names(de_list)){
  res <- de_list[[strat]]
  res[row.names(names_key), 'gene_name_plot'] <- names_key$gene_name_plot
  #selectLab <- res[which(res$padj <= pthresh & abs(res$log2FoldChange) >= FCthresh), 'gene_name_plot']
  selectLab <- res[which(res$padj <= 0.01 & abs(res$log2FoldChange) >= 1.5 & !grepl('^ENSTGUG', res$gene_name_plot)), 'gene_name_plot']
  #selectLab <- res[which(res$padj <= 0.01 & abs(res$log2FoldChange) >= 3), 'gene_name_plot']
  #selectLab <- res[which(res$padj <= pthresh & abs(res$log2FoldChange) >= FCthresh & !grepl('^ENSTGUG', res$gene_name_plot)), 'gene_name_plot']
  
  p3 <- EnhancedVolcano(res,
                       lab = res[,'gene_name_plot'],
                       selectLab = selectLab,
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = strat,
                       caption = '',
                       subtitle = paste0(test_desc, ' DE Volcano Plot'),
                       xlab = bquote(~Log[2]~ 'fold change'),
                       pCutoff = pthresh,
                       FCcutoff = FCthresh,
                       pointSize = 2.0,
                       labSize = 2.5,
                       labhjust = 0.5,
                       labvjust = 1.5,
                       colAlpha = 1,
                       legendPosition = 'top',
                       #boxedLabels = T,
                       drawConnectors = TRUE,
                       widthConnectors = 0.5,
                       lengthConnectors = unit(0.02, "npc"),
                       arrowheads = F
  )
  p3 <- p3 + theme(legend.position="none")
  print(p3)
}

### Manually add some 
additional_genes <- c('SLC26A11', 'ROS1', 'ENKUR')
names(additional_genes) <- c('ENSTGUG00000003079', 'ENSTGUG00000011909', 'ENSTGUG00000001164')

genes <- c(genes, additional_genes)
### Boxplots of normalized counts, with DE pvalue and logFC in title
for(id in names(genes)){
  name <-  gtf[id, 'gene_name']
  plot_title <- ifelse(is.na(name), id, name)
  for(set in names(de_list)){
    if(set != 'All'){
      plot_dat <- as.data.frame(covs[which(covs[, strat_title] == set),])
    } else{
      plot_dat <- as.data.frame(covs)
    }
    plot_dat[, 'norm_counts'] <- as.numeric(norm_counts[id, row.names(plot_dat)])
    adjp <- signif(de_list[[set]][id, 'padj'], digits = 3)
    l2fc <- signif(de_list[[set]][id, 'log2FoldChange'], digits = 3)
    p4 <- ggplot(plot_dat, aes(x = get(test), y = norm_counts)) + 
      labs(title = paste0(plot_title, ' (', set, ')'),
           subtitle = paste0(test_desc, ' DESeq2 padj: ', adjp, ', log2FC: ', l2fc)) + xlab(names(test)) +
      geom_boxplot()
    ### Show all points if there are few enough
    if(dim(plot_dat)[1] < 50){
      p4 <- p4 + geom_jitter(color = "red", size = 1, alpha = 0.9)
    }
    print(p4)
  }
}

dev.off()
