#library('sys')
#library('viridis')
library('data.table')
#library('PKI')
library('stringr')
library('stringi')

source('/home/cwake/snakemakes/sc_functions.R')
source('/home/cwake/snakemakes/Utility_functions.R')

if(interactive()){
  gtf_file <- '/home/cwake/resources/genomes/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
  rds_file <- '/home/cwake/data/gtf.RDS'
} else{
  args = commandArgs(trailingOnly=TRUE)
  ### snakemake input
  gtf_file <- args[1]
  rds_file <- args[2]
}

gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
saveRDS(gtf, rds_file)
