#library('sys')
library('DESeq2')
library('ggplot2')
#library('stringr')
#library('GenomicRanges')
library('dplyr')

source('/home/cwake/snakemakes/Utility_functions.R')
#source('/home/cwake/snakemakes/DE_functions.R')

args = commandArgs(trailingOnly=TRUE)
count_in <- args[1]
covs_in <- args[2]
candidates <- args[3]
method <- args[4]
txt_out <- args[5]
pdf_out <- args[6]

# count_in <- '/home/cwake/projects/2021612_finch/data/counts/normalizedCounts.txt'
# covs_in <- '/home/cwake/projects/2021612_finch/data/Covariates_QC_metrics_filter.csv'
# candidates <- 'plate, RNA_Extraction_Batch'
# method <- 'model'
# txt_out <- '/home/cwake/projects/2021612_finch/data/batch_evaluation.txt'
# pdf_out <- '/home/cwake/projects/2021612_finch/data/batch_evaluation.pdf'

### Analysis only done in the 'else' statement, if the user has input some candidate batch variable to evaluate. If they haven't, only output a message recommending that they do.
if(candidates == ''){
  cat(paste0('Automated recommendation: Enter batch information in config file', ''), file = txt_out, sep = "\n", append = T)
  cat('Model batch: ', file = txt_out, sep = "\n", append = T)
  cat(paste0('Automated recommendation: Enter batch information in config file', ''), file = pdf_out, sep = "\n")
} else {
  ### Set default method
  if(method == ''){
    method <- 'model'
  }
  candidates <- trimws(strsplit(candidates,',')[[1]])
  ### Read counts
  norm_counts <- read.csv(count_in, stringsAsFactors=F, header = T,
                          check.names = F, 
                          colClasses = 'character',
                          sep = '\t')
  ### Convert to numeric
  norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))
  
  ### Read in covariates
  covs <- read.csv(covs_in, stringsAsFactors = F, header = T,
                   check.names = F, 
                   colClasses = 'character',
                   sep = ',')
  covs <- covs[which(covs$filter == '0'),]
  
  PCA <- run_pca(norm_counts)
  pcs <- c(1,2)
  note <- ''
  pc1 <- paste0('PCA',pcs[1])
  pc2 <- paste0('PCA',pcs[2])
  pca <- as.data.frame(PCA$x)
  ### Reorder phase to match pca
  phases <- covs[match(rownames(pca), covs$ID), ]
  
  covs_numeric <- c('N_counts', 'Aligned_input', 'Aligned_uniquely', 'Aligned_multi', 'RIN_score')
  covs_factor <- c('Euth_Age','Sex','Treatment','Cohort')
  ### Convert to numeric and factor
  covs[, covs_numeric] <- mutate_all(covs[, covs_numeric], function(x) as.numeric(x))
  covs[, covs_factor] <- mutate_all(covs[, covs_factor], function(x) as.factor(x))
  
  pdf(pdf_out)
  for(b in candidates){
    lvls <- unique(covs[,b])
    ### Separate numeric from non-numeric, to sort the latter
    lvls <- c(lvls[which(is.na(suppressWarnings(as.numeric(lvls))))], as.character(sort(as.numeric(lvls[which(!is.na(suppressWarnings(as.numeric(lvls))))]))))
    covs[, b] <- factor(covs[, b], levels = lvls)
    table(covs[, b])
    ###### BASICS
    ### PCA plot
    df <- data.frame(x = pca[,pcs[1]], y = pca[,pcs[2]],
                     phase = phases[, b])
    
    p2 <- ggplot(df, aes(x = x, y = y, colour = phase)) +
      geom_point() +
      theme_bw() +
      xlab(pc1) + ylab(pc2) + ggtitle(paste0(pc1, ' vs. ', pc2, ' ', note))
    p2
    
    ### Boxplot of per-sample counts
    boxplot(formula = log(Aligned_uniquely) ~ eval(parse(text = b)),
            data = covs, main = paste0(b,' Batches log aligned N'), 
            xlab = 'Batches',ylab = 'log N reads',col='green',outcol ='red')
    batches = levels(covs[, b])
    
    ### Customize axes start and text size based on the number of batches (assuming 4 to 11). Might need to update with more batches.
    text_x = 0 + (11 - length(batches))*(.4/7)
    text_size = 0.75 + (11 - length(batches))*(.25/7)
    text(x = text_x, y = 17.0, 'Mean',pos = 4, cex = text_size)
    text(x = text_x, y = 16.9, 'Median', pos = 4, cex = text_size)
    text(x = text_x, y = 16.8, 'sd', pos = 4, cex = text_size)
    for(i in batches){
      batch_mean = as.character(round(mean(log(covs[which(covs[, b] == i),]$Aligned_uniquely)), digits = 2))
      batch_median = as.character(round(median(log(covs[which(covs[, b] == i),]$Aligned_uniquely)), digits = 2))
      batch_sd = as.character(round(sd(log(covs[which(covs[, b] == i),]$Aligned_uniquely)), digits = 2))
      text(x = i, y = 17, batch_mean, pos = 4, cex = text_size)
      text(x = i, y = 16.9, batch_median, pos = 4, cex = text_size)
      text(x = i, y = 16.8, batch_sd, pos = 4, cex = text_size)
    }
    
    ### Test batch for association with covs_numeric and covs_factor, returning p-values   
    cov_assoc <- batch_cov_associations(covs, b, covs_numeric, covs_factor)
    
    #
    
  }
  #pdf(pdf_out)
  out <- gene_cov_associations(norm_counts, covs, covs_batch = candidates)
  res_p <- out[['p']]
  res_c <- out[['t']]
  ### Create FDR adjusted p values
  res_padj <- res_p
  for(cov in colnames(res_p)){
    res_padj[,cov] <- p.adjust(res_p[,cov], method = "fdr")
  }
  #pdf(pdf_out)
  for(i in 1:length(colnames(res_p))){
    p1 <- ggplot(res_padj, aes(x = eval(parse(text = colnames(res_p)[i])))) + 
      geom_histogram() + ggtitle(paste0('adjusted p-values, ', colnames(res_p)[i])) + xlab('p') + ylab('feature count')
    p2 <- ggplot(res_c, aes(x = eval(parse(text = colnames(res_p)[i])))) + 
      geom_histogram() + ggtitle(paste0('statistic, ', colnames(res_p)[i])) + xlab('stat') + ylab('feature count')
    print(p1)
    print(p2)
  }
  dev.off()
  
  txt_summary <- genewise_summarize(res_p, res_padj, res_c)
  write.table(txt_summary, txt_out, quote = F, sep = '\t', row.names = T, col.names = T)
  
  if(any(txt_summary$`Fraction p < 0.05` > 0.1)){
    bestbatch <- row.names(txt_summary[which(txt_summary$`Fraction p < 0.05` == max(txt_summary$`Fraction p < 0.05`)),])
  } else{
    bestbatch <- ''
  }
  cat(paste0('Automated recommendation for model: ', bestbatch), file = txt_out, sep = "\n", append = T)
  cat(paste0('Automated recommendation for ComBat: ', ''), file = txt_out, sep = "\n", append = T)
  cat(paste0('Model batch: ', bestbatch), file = txt_out, sep = "\n", append = T)
  cat(paste0('ComBat batch: ', ''), file = txt_out, sep = "\n", append = T)
}

