library(spacexr)
library(Matrix)
library(doParallel)

source('../iter_optim.R')
source('../initialization.R')
source('../clustering.R')
source('../analysis.R')

puck <- readRDS('../puckCropped.rds')
cell_types <- c(1:10)

myRCTD <- create.RCTD.noref(puck, cell_type_names = cell_types)
myRCTD <- set_internal_vars(myRCTD)
myRCTD <- assign.cell.types(myRCTD, runif_doublet)
myRCTD@config$gene_cutoff <- 0; myRCTD@config$fc_cutoff <- 0; myRCTD@config$gene_cutoff_reg <- 0; myRCTD@config$fc_cutoff_reg <- 0
#myRCTD@config$max_cores <- 1
RCTD_results <- run.iter.optim(myRCTD, cell_types = cell_types, n_iter = 5)


saveRDS(RCTD_results,'../RCTDlist_randominit.rds')
RCTD_results <- readRDS('../RCTDlist_randominit.rds')
