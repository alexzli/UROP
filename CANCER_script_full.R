datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c( .libPaths(), file.path(datadir, 'Rpcg')))
library(spacexr)
library(Matrix)
source(file.path(datadir, '/R/optimization.R'))
source(file.path(datadir, '/R/initialization.R'))

# data loading
myRCTD <- readRDS(file.path(datadir, '/data/myRCTD_201014_03.rds'))
results <- read.csv(file = file.path(datadir, '/data/results_processed_slideseq_data_2020-12-12_Puck_201014_03_2020-12-12_Puck_201014_03_cell_types.csv'))
puck <- readRDS(file.path(datadir, '/data/cancer_puck.rds'))

# data preprocessing
original_info <- get_norm_ref(puck, myRCTD@cell_type_info$info[[1]], rownames(myRCTD@cell_type_info$info[[1]]), myRCTD@internal_vars$proportions)  
cancer_results <- results[results$HRS_location == TRUE,]
cell_types <- rep('ignore', length(colnames(puck@counts)))
names(cell_types) <- colnames(puck@counts)
for (barcode in names(cell_types)) {
  if (barcode %in% cancer_results$barcodes)
    cell_types[barcode] = 'Cancer_cells'
}
cell_types <- factor(cell_types)
cancer_info <- get_cell_type_info(puck@counts, cell_types, puck@nUMI)[[1]]

genes <- intersect(rownames(cancer_info),rownames(original_info))
cancer_info <- cancer_info[1:1]
cancer_info <- cancer_info[genes,]
original_info <- original_info[genes,]
original_info$Cancer_cells <- cancer_info
original_info[is.na(original_info)] <- 0
#original_info <- sweep(original_info, 2, colSums(original_info), '/')
new_info <- list(original_info, append(myRCTD@cell_type_info$info[[2]], 'Cancer_cells'), myRCTD@cell_type_info$info[[3]]+1)
cell_type_info <- list(info = new_info, renorm = new_info)
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
myRCTD@config$max_cores <- 32

RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, CELL_MIN_INSTANCE = 0)
saveRDS(RCTDlist, file.path(datadir, '/data/RCTD_list_cancer_FULL.rds'))