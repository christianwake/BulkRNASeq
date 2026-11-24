
#library('sys')
library('dplyr')
#library('viridis')
library('readxl')
library('stringr')
library('WriteXLS')
library('data.table')
library('tools')
source('/home/cwake/snakemakes/sample_sheet_functions.R')
#source('/home/cwake/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021612_finch'
  investigator <- 'Rachel_Davis'
  reference <- 'tguttata'
  run_path <- '/home/cwake/data/'
  csv_file <- paste0('/home/cwake/projects/', project, 
                     '/SampleSheets/Sample_sheet.csv')
  ## bcl2fastq sample sheets from excel docs
  #excel_file <- '/home/cwake/projects/2021612_finch/RD_TLbulkrnaseq_052721.xlsx'
}else{
  args = commandArgs(trailingOnly=TRUE)
  project <- args[1]
  investigator <- args[2]
  reference <- args[3]
  run_path <- args[4]
  csv_file <- args[5]
}

### Read csv
dat <- read.csv(csv_file, quote = "", header = T, check.names = F)
colnames(dat)[which(colnames(dat) == 'index')] <- 'RC'

### Conforming column names and other formatting
dat <- sample_file_to_csv(dat, bcl2fastq2 = F)

### Full flowcell paths using ls and the abbreviated name
if(!'flowcell_full' %in% colnames(dat)){
  full_flowcell <- sapply(unique(dat$FCID), function(x) list.files(path = paste0(run_path), pattern = paste0('*', x), include.dirs = T)[1])
  dat$flowcell_full <- full_flowcell[dat$FCID]
} else{
  full_flowcell <- unique(dat$flowcell_full)
  names(full_flowcell) <- sapply(full_flowcell, function(x) unique(dat[which(dat$flowcell_full == x), 'FCID']))
}
if(!'Sample_Project' %in% colnames(dat)){
  dat[, 'Sample_Project'] <- project
}
##### Now per flowcell bcl2fastq format csv
### Get unique flowcells from dat
a <- unique(dat[, c('FCID', 'flowcell_full')])
flowcells <- a$flowcell_full
names(flowcells) <- a$FCID

### bcl2fastq format header
header <- paste0('[Header]\nDate,', Sys.Date(), '\nWorkflow,bcl2fastq2\nInvestigator,', investigator, '\n[Data]')
#### For each flowcell, write an individual csv file
for(i in 1:length(flowcells)){
  ### Subset rows
  fc_dat <- dat[which(dat$FCID == names(full_flowcell)[i]), ]
  ### Subset columns and change format
  fc_dat <- bcl2fastq2_format(fc_dat, indeces = c('RC', 'i7_index'), Lane = 'Lane')
  file_path <- paste0(run_path, full_flowcell[i], '/', 'SampleSheet.csv')
  write(header, file = file_path, append = F)
  suppressWarnings(write.table(fc_dat, file_path, quote = F, sep = ',', col.names = T, row.names = F, append = T))
}

