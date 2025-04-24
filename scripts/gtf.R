### This script assumes you've already completed bcl2fastq, Trimmomatic, STAR and maybe BALDR
library('sys')
library('viridis')
library('data.table')
library('PKI')
library('stringr')
library('stringi')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  gtf_file <- '/data/vrc_his/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
  rds_file <- '/data/vrc_his/douek_lab/wakecg/data/gtf.RDS'
} else{
  args = commandArgs(trailingOnly=TRUE)
  ### snakemake input
  gtf_file <- args[1]
  rds_file <- args[2]
}

gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
saveRDS(gtf, rds_file)
