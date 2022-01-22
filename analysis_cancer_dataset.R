source('../UROP/R/optimization.R')
source('../UROP/R/initialization.R')
source('../UROP/R/clustering.R')
source('../UROP/R/analysis.R')

myRCTD <- readRDS('../UROP/data/myRCTD_201014_03.rds')
results <- read.csv(file = '../UROP/data/results_processed_slideseq_data_2020-12-12_Puck_201014_03_2020-12-12_Puck_201014_03_cell_types.csv')

puck <- myRCTD@originalSpatialRNA
original_info <- myRCTD@cell_type_info$info
cancer_results <- results[results$HRS_location == TRUE,]
cell_types <- rep('ignore', length(colnames(puck@counts)))
names(cell_types) <- colnames(puck@counts)
for (barcode in names(cell_types)) {
  if (barcode %in% cancer_results$barcodes)
    cell_types[barcode] = 'Cancer_cells'
}
cancer_info <- get_cell_type_info(puck@counts, factor(cell_types), puck@nUMI)[[1]]

new_info <- as.data.frame(matrix(0, nrow=length(rownames(original_info[[1]])), ncol=1))
rownames(new_info) <- rownames(original_info[[1]])
colnames(new_info) <- 'Cancer_cells'
for (gene in rownames(original_info[[1]])) {
  if (gene %in% rownames(cancer_info))
    new_info[gene,'Cancer_cells'] = cancer_info[gene,'Cancer_cells']
}

new_info <- list(cbind(original_info[[1]], new_info), append(original_info[[2]], 'Cancer_cells'), original_info[[3]]+1)
cell_type_info <- list(info = new_info, renorm = NULL)
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, CELL_MIN_INSTANCE = 10, used_reference = TRUE)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_cancer.rds')

