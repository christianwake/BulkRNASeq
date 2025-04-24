library('readxl')

file <- 'C:/Users/wakecg/Documents/SCSeq/ZebraFinch/DNA_RNA Extractions.xlsx'
dat <- as.data.frame(read_xlsx(file, sheet = 'Bone Marrow RNA'))

source('C:/Users/wakecg/Documents/SCSeq/sc_functions.R')

plots_path <- 'C:/Users/wakecg/Documents/SCSeq/ZebraFinch/plots/'
covs_numeric <- c()
covs_factor <- c('Individual', 'Treatment', 'Cohort', 'Euth Assignment', 'Date of sample',
                  'Extraction Batch')
cov_associations <- pairwise_associations(dat, pc = 'p', plots_path = plots_path,
                                           covs_numeric = covs_numeric, 
                                           covs_factor = covs_factor,
                                           covs_binary = c())
adjps <- cov_associations[['res']]
test_names <- cov_associations[['test']]
