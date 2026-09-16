library('sys')
library('readr')
library('WriteXLS')
library('data.table')

source('/home/cwake/snakemakes/Utility_functions.R')
source('/home/cwake/snakemakes/DE_functions.R')

get_strats <- function(gs, dats, n = F){
  strats <- c()
  ### For each stratification
  for(i in 1:length(dats)){
    if(gs %in% row.names(dats[[i]])){
      strats <- c(strats, names(dats)[i])
    }
  }
  if(n){
    res <- length(strats)
  } else{
    res <- toString(strats)
  }
  return(res)
}

if(interactive()){
  args <- c('/home/cwake/projects/2021612_finch/results/Treatment_DESeq2.xls',
            '/home/cwake/projects/2021612_finch/results/Treatment_DESeq2_sig.xls',
            '0.05',
            '/home/cwake/projects/2021612_finch/results/Treatment/All/DESeq2_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D6/DESeq2_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D18/DESeq2_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D35/DESeq2_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D90/DESeq2_results.txt')
  
  args <- c('/home/cwake/projects/2021612_finch/results/Treatment_GSEA.xls',
            '/home/cwake/projects/2021612_finch/results/Treatment_GSEA_sig.xls',
            '0.05',
            '/home/cwake/projects/2021612_finch/results/Treatment/All/fgsea_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D6/fgsea_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D18/fgsea_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D35/fgsea_results.txt',
            '/home/cwake/projects/2021612_finch/results/Treatment/Euth_Age-D90/fgsea_results.txt')

}else{
  args = commandArgs(trailingOnly=TRUE)
}
xls_file <- args[1]
sig_file <- args[2]
pthresh <- as.numeric(args[3])
files <- args[4:length(args)]

dat_list <- vector(mode = "list", length = length(files))
tab_names <- c()
for(i in 1:length(files)){
  file <- files[i]
  ### Get test and strat names from file name
  a <- tail(strsplit(dirname(file), '/')[[1]], 2)
  tab_names <- c(tab_names, paste0(a[1], '_', a[2]))
  dat_list[[i]] <- read.table(file, sep = '\t', quote = "", header = T, row.names = 1)
}
names(dat_list) <- tab_names
dats <- lapply(dat_list, function(dat) dat[which(dat$padj < as.numeric(pthresh)),])

### Create summary data
sig_res <- unique(unlist(lapply(dats, function(dat) row.names(dat))))
out <- as.data.frame(sig_res)
row.names(out) <- sig_res
### Add size, gtf_name
cols <- c('size', 'gene_name_gtf')
cols <- cols[cols %in% colnames(dat_list[[1]])]
for(col in cols){
  out[, col] <- sapply(row.names(out), function(x) dat_list[[1]][x, col])
}
out$n_sig <- sapply(row.names(out), function(gs) get_strats(gs, dats, T))
out$sig_strats <- sapply(row.names(out), function(gs) get_strats(gs, dats))
### Specific to finch data
# out$sig_strats <- gsub('Treatment_Euth_Age-', '', out$sig_strats)
# out$sig_strats <- gsub('Treatment_', '', out$sig_strats)

out <- out[order(out$n_sig, decreasing = T),]
out <- out[, colnames(out)[colnames(out) != 'sig_res']]

dat_list[['Summary']] <- out
dats[['Summary']] <- out
### Write all results
WriteXLS(x = dat_list, ExcelFileName = xls_file, SheetNames = names(dat_list), row.names = T)
### Write significant results
WriteXLS(x = dats, ExcelFileName = sig_file, SheetNames = names(dats), row.names = T)
###

