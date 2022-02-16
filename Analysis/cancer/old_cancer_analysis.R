source('../UROP/R/optimization.R')
source('../UROP/R/initialization.R')
source('../UROP/R/clustering.R')
source('../UROP/R/analysis.R')

# Initial code for analyzing cancer dataset. Updated code can be found in run_cancer.R,
# and updated analysis can eb found in analyze_cancer.R


## --------------------------------- DATA PREPROCESSING ---------------------------------

# load initial results and assignments, cancer puck
myRCTD <- readRDS('../UROP/data/myRCTD_201014_03.rds')
results <- read.csv(file = '../UROP/data/results_processed_slideseq_data_2020-12-12_Puck_201014_03_2020-12-12_Puck_201014_03_cell_types.csv')
puck <- readRDS('../UROP/data/cancer_puck.rds')


# load full gene list (not currently functional)
#counts <- read.csv(file = '../UROP/data/MappedDGEForR.csv')
#rownames(counts) <- counts$GENE; counts$GENE <- NULL
#nUMI <- colSums(counts)
#coords <- myRCTD@originalSpatialRNA@coords
#puck <- SpatialRNA(coords, counts, nUMI)
#saveRDS(puck, '../UROP/data/cancer_puck.rds')


# extract cancer cell type info
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

## --------------------------------- INITIAL DATA ANALYSIS ---------------------------------

# spatially plot initial cancer cell types
results <- myRCTD@results$results_df
results$first_type = factor(results$first_type, levels=append(levels(results$first_type), 'Cancer_cells'))
results$second_type = factor(results$second_type, levels=append(levels(results$second_type), 'Cancer_cells'))
for (barcode in cancer_results$barcodes) {
  results[barcode,'first_type'] = 'Cancer_cells'
  results[barcode,'second_type'] = 'Cancer_cells'
}
plot_all_cell_types(results, myRCTD@originalSpatialRNA@coords, levels(results$first_type), '..')


# gene expression dimensionality reduction plots (UMAP/tSNE/PCA)
slide.seq <- CreateSeuratObject(counts = puck@counts)
slide.seq <- NormalizeData(slide.seq)
slide.seq <- FindVariableFeatures(slide.seq)
slide.seq <- ScaleData(slide.seq)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- RunTSNE(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq.clusters <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
for (barcode in names(slide.seq.clusters@active.ident)) {
  if (barcode %in% cancer_results$barcodes) {
    slide.seq.clusters@active.ident[barcode] = 1
  }
  else {
    slide.seq.clusters@active.ident[barcode] = 0
  }
}
plot1 <- DimPlot(slide.seq.clusters, reduction = "tsne")
plot1


## --------------------------------- SEMISUPERVISED CANCER CELL TYPE PREDICTION ALGORITHM ---------------------------------
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
class_df <- myRCTD@internal_vars$class_df
class_df['Cancer_cells', 'class'] = 'Cancer_cells'
myRCTD@internal_vars$class_df <- class_df
myRCTD@config$max_cores <- 4

myRCTD <- fit.gene.expression(myRCTD, cell_types = cell_types, CELL_MIN_INSTANCE = 0, sigma_gene = FALSE)
#RCTD <- fitPixels(myRCTD, doublet_mode = 'doublet')

RCTDlist <- run.iter.optim(myRCTD, cell_types = cell_types)
saveRDS(RCTDlist, '../UROP/objects/RCTD_list_cancer_de_info.rds')
        


## --------------------------------- IGNORE BELOW ---------------------------------

# old semisupervised algorithm using mean info initially
genes <- intersect(rownames(cancer_info),rownames(original_info))
cancer_info <- cancer_info[1:1]
cancer_info <- cancer_info[genes,]
original_info <- original_info[genes,]
original_info$Cancer_cells <- cancer_info
original_info[is.na(original_info)] <- 0
#original_info <- sweep(original_info, 2, colSums(original_info), '/')
new_info <- list(original_info, append(myRCTD@cell_type_info$info[[2]], 'Cancer_cells'), myRCTD@cell_type_info$info[[3]]+1)
cell_type_info <- list(info = new_info, renorm = new_info)
#gene_list_reg <- myRCTD@internal_vars$gene_list_reg
#gene_list_bulk <- myRCTD@internal_vars$gene_list_bulk
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
#myRCTD@config$RCTDmode <- "doublet"
#myRCTD <- choose_sigma_c(myRCTD)
#myRCTD <- fitPixels(myRCTD, doublet_mode = "doublet")
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, CELL_MIN_INSTANCE = 0)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_cancer_FULL.rds')

# initial results analysis
RCTDlist <- readRDS('../UROP/objects/RCTD_list_cancer.rds')
RCTDlist <- readRDS('../UROP/objects/RCTD_list_cancer_FULL.rds')
RCTDlist <- readRDS('../UROP/objects/RCTD_list_cancer_de_info.rds')
RCTDinit <- readRDS('../UROP/objects/RCTD_init_cancer.rds')
for (i in 1:length(RCTDlist)) {
  print(cell.confusion.mat(RCTDinit, RCTDlist[[i]])$overall['Accuracy'])
}
for (i in 1:(length(RCTDlist)-1)) {
  print(cell.confusion.mat(RCTDlist[[i]], RCTDlist[[i+1]])$overall['Accuracy'])
}
for (i in 1:(length(RCTDlist)-1)) {
  print(gene.mse(RCTDlist[[i]], RCTDlist[[i+1]]))
}
head(RCTDlist[[8]]@cell_type_info$renorm[[1]][order(RCTDlist[[8]]@cell_type_info$renorm[[1]]$Cancer_cells, decreasing=TRUE),]['Cancer_cells'], 10)
table(RCTDlist[[8]]@results$results_df[RCTDlist[[8]]@results$results_df$spot_class == 'singlet',]$first_type)
plot.cell.types(RCTDlist[[8]])

# binary cancer classification
new_info <- list(cancer_info, c('ignore', 'Cancer_cells'), 2)
cell_type_info <- list(info = new_info, renorm = NULL)
gene_list_reg <- myRCTD@internal_vars$gene_list_reg
gene_list_bulk <- myRCTD@internal_vars$gene_list_bulk
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info, gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk)
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = TRUE)
saveRDS(RCTDlist, '../UROP/objects/RCTD_list_cancer_binary.rds')
RCTDlist <- readRDS('../UROP/objects/RCTD_list_cancer_binary.rds')
# multiclass cancer classification
RCTD <- readRDS('../UROP/data/myRCTD_201014_03.rds')
doublet_weights <- RCTD@results$weights_doublet
cancer_results <- results[results$HRS_location == TRUE,]
cancer_weights <- doublet_weights[intersect(rownames(doublet_weights),cancer_results$barcodes),]
rownames(cancer_results) <- cancer_results$barcodes
cancer_results <- cancer_results[,c('spot_class', 'first_type', 'second_type')]
cancer_results$spot_class <- factor(cancer_results$spot_class, levels=c("reject", "singlet", "doublet_certain", "doublet_uncertain"))
cancer_results$first_type <- factor(cancer_results$first_type, levels=cell_types)
cancer_results$second_type <- factor(cancer_results$second_type, levels=cell_types)
myRCTD <- create.RCTD.noref(puck)
myRCTD <- set_internal_vars(myRCTD)
myRCTD@results$results_df <- cancer_results
myRCTD@results$weights_doublet <- cancer_weights
myRCTD@internal_vars$cell_types_assigned <- TRUE
myRCTD@config$RCTDmode <- "doublet"
myRCTD@cell_type_info <- list(info = list(data.frame(), cell_types, length(cell_types)), renorm = NULL)
myRCTD <- fit.gene.expression(myRCTD, cell_types, CELL_MIN_INSTANCE = 0)
saveRDS(myRCTD,'../UROP/objects/RCTD_cancer_expression.rds')
