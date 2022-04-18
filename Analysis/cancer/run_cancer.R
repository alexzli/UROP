datadir <- '~/UROP/'
source(file.path(datadir, 'R/algorithm.R'))


# LOAD DATA
myRCTD <- readRDS(file.path(datadir, '/Data/myRCTD_201014_03.rds'))
results <- read.csv(file = file.path(datadir, '/Data/results_processed_slideseq_data_2020-12-12_Puck_201014_03_2020-12-12_Puck_201014_03_cell_types.csv'))


# DATA PREPROCESSING
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
class_df <- myRCTD@internal_vars$class_df
class_df['Cancer_cells', 'class'] = 'Cancer_cells'
myRCTD@internal_vars$class_df <- class_df
myRCTD@config$max_cores <- 4


# FULL GENE LIST
counts <- read.csv(file = file.path(datadir, '/Data/MappedDGEForR.csv'))
rownames(counts) <- counts$GENE; counts$GENE <- NULL
nUMI <- colSums(counts)
myRCTD <- readRDS(file.path(datadir, '/Data/myRCTD_201014_03.rds'))
coords <- myRCTD@originalSpatialRNA@coords
puck <- SpatialRNA(coords, counts, nUMI)
saveRDS(puck, file.path(datadir, '/Data/cancer_puck.rds'))


saveRDS(myRCTD, file.path(datadir, '/Objects/cancer_RCTD_initial.rds'))


# RUN ALGORITHM
myRCTD <- readRDS(file.path(datadir, 'Objects/cancer_RCTD_initial.rds'))
RCTD_list <- run.semisupervised(myRCTD)
saveRDS(RCTD_list, file.path(datadir, '/Objects/cancer_RCTD_final.rds'))
