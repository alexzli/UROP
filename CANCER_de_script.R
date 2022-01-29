datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c( .libPaths(), file.path(datadir, 'Rpcg')))
library(spacexr)
library(Matrix)
source(file.path(datadir, '/R/optimization.R'))
source(file.path(datadir, '/R/initialization.R'))

myRCTD <- readRDS(file.path(datadir, '/data/myRCTD_201014_03.rds'))
results <- read.csv(file = file.path(datadir, '/data/results_processed_slideseq_data_2020-12-12_Puck_201014_03_2020-12-12_Puck_201014_03_cell_types.csv'))

cancer_results <- results[results$HRS_location == TRUE,]
results <- myRCTD@results$results_df
results$first_type = factor(results$first_type, levels=append(levels(results$first_type), 'Cancer_cells'))
results$second_type = factor(results$second_type, levels=append(levels(results$second_type), 'Cancer_cells'))
for (barcode in cancer_results$barcodes) {
  results[barcode, 'spot_class'] = 'singlet'
  results[barcode,'first_type'] = 'Cancer_cells'
}
myRCTD@results$results_df <- results
myRCTD@cell_type_info$info[[2]] <- append(myRCTD@cell_type_info$info[[2]], 'Cancer_cells')
myRCTD@cell_type_info$renorm[[2]] <- append(myRCTD@cell_type_info$renorm[[2]], 'Cancer_cells')
myRCTD@cell_type_info$info[[3]] <- myRCTD@cell_type_info$info[[3]] + 1
myRCTD@cell_type_info$renorm[[3]] <- myRCTD@cell_type_info$renorm[[3]] + 1
cell_types <- levels(results$first_type)
myRCTD@config$max_cores <- 32

RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = TRUE, cell_types = cell_types)
saveRDS(RCTDlist, file.path(datadir, '/data/RCTD_list_cancer_degeneratedinfo.rds')
