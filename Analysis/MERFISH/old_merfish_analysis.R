library(STdeconvolve)
source('../UROP/R/ITERATIVE_OPTIMIZATION.R')
source('../UROP/R/analysis.R')

# Initial code for analyzing old MERFISH dataset. Updated code can be found in
# analyze_merfish.R.


puck <- readRDS('../UROP/objects/MERFISH_puck_20um2.rds')

# PLOT EMBEDDINGS

slide.seq <- CreateSeuratObject(counts = puck@counts)
slide.seq <- NormalizeData(slide.seq)
slide.seq <- FindVariableFeatures(slide.seq)
slide.seq <- ScaleData(slide.seq)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- RunTSNE(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq.clusters <- FindClusters(slide.seq, resolution = 0.18, verbose = FALSE)
slide.seq.clusters@active.ident <- factor(unlist(true_types))
slide.seq.clusters@active.ident <- pred_types
plot1 <- DimPlot(slide.seq.clusters, reduction = "tsne", label = TRUE)
plot1

pred_types <- RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet',]$first_type
names(pred_types) <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet',])
true_types <- as.data.frame(simFN7$gtSpotTopics)
true_types <- apply(true_types, 1, function(row_tmp){
  names(true_types)[which(row_tmp == 1)]
})

# IGNORE
assignments <- gen.clusters(puck, resolution = 0.15, SCT = F)
myRCTD <- create.RCTD.noref(puck)
myRCTD <- set_internal_vars(myRCTD)
myRCTD <- assign.cell.types(myRCTD, assignments)

true_barcodes <- rownames(simFN7$cellCounts[simFN7$cellCounts$counts == 1,])
pred_barcodes <- rownames(myRCTD@results$results_df[myRCTD@results$results_df$spot_class == 'singlet',])
barcodes <- intersect(true_barcodes, pred_barcodes)
pred_types <- myRCTD@results$results_df[barcodes,'first_type']
true_types <- as.data.frame(simFN7$gtSpotTopics[barcodes,])
true_types <- apply(true_types, 1, function(row_tmp){
  names(true_types)[which(row_tmp == 1)]
})
mytable <- table(true_types, pred_types)
assignment.accuracy(mytable)
assignment.accuracy.plot(mytable)
mytable


## ANALYZE RESULTS
RCTD_list <- readRDS('../UROP/results/MERFISH_results/MERFISH20um2_results_mean_fit_2_10.rds')
RCTDpred <- RCTD_list[[length(RCTD_list)]]
simFN7 <- readRDS('../UROP/objects/simFN7_20um2.rds')

plot.cell.types(RCTDpred)
table(RCTDpred@results$results_df$spot_class)
table(simFN7$cellCounts$counts)

rejected_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'reject',])
table(simFN7$cellCounts[rejected_barcodes,]$counts)
singlet_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet',])
table(simFN7$cellCounts[singlet_barcodes,]$counts)
doublet_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'doublet_certain',])
table(simFN7$cellCounts[doublet_barcodes,]$counts)
true_doublet_barcodes <- rownames(simFN7$cellCounts[simFN7$cellCounts$counts == 2,])
table(RCTDpred@results$results_df[true_doublet_barcodes,]$spot_class)

# table of number of cells
common_barcodes = intersect(rownames(RCTDpred@results$results_df), rownames(simFN7$cellCounts))
table(RCTDpred@results$results_df[common_barcodes,'spot_class'], simFN7$cellCounts[common_barcodes, 'counts'])

num_correct_singlets = length(intersect(rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet',]), rownames(simFN7$cellCounts[simFN7$cellCounts$counts == 1,])))
num_correct_conf_doublets = length(intersect(rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'doublet_certain',]), rownames(simFN7$cellCounts[simFN7$cellCounts$counts == 2,])))
(num_correct_singlets + num_correct_conf_doublets) / length(rownames(RCTDpred@results$results_df))
num_correct_singlets/length(rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet',]))
num_correct_conf_doublets/length(rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'doublet_certain',]))


# SINGLET ANALYSIS
true_barcodes <- rownames(simFN7$cellCounts[simFN7$cellCounts$counts == 1,])
pred_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet',])
barcodes <- intersect(true_barcodes, pred_barcodes)

pred_types <- RCTDpred@results$results_df[barcodes,'first_type']
true_types <- as.data.frame(simFN7$gtSpotTopics[barcodes,])
true_types <- apply(true_types, 1, function(row_tmp){
  names(true_types)[which(row_tmp == 1)]
})
mytable <- table(true_types, pred_types)
assignment.accuracy(mytable)
assignment.accuracy.plot(mytable)
mytable

for (i in 1:length(RCTD_list)) {
  RCTDpred <- RCTD_list[[i]]
  true_barcodes <- rownames(simFN7$cellCounts[simFN7$cellCounts$counts == 1,])
  pred_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet',])
  barcodes <- intersect(true_barcodes, pred_barcodes)
  
  pred_types <- RCTDpred@results$results_df[barcodes,'first_type']
  true_types <- as.data.frame(simFN7$gtSpotTopics[barcodes,])
  true_types <- apply(true_types, 1, function(row_tmp){
    names(true_types)[which(row_tmp == 1)]
  })
  mytable <- table(true_types, pred_types)
  print(mytable[c('Excitatory', 'Inhibitory'),c('0','7')])
}

# GENE EXPRESSION ANALYSIS
head(simFN7$gtCtGenes)

ref_expression <- data.matrix(log(t(simFN7$gtCtGenes)))
pred_expression <- data.matrix(log(RCTDpred@cell_type_info$renorm[[1]]))
gene_list <- intersect(rownames(ref_expression), rownames(pred_expression))
ref_expression <- ref_expression[gene_list,]
pred_expression <- pred_expression[gene_list,]
ref_cells <- colnames(ref_expression); pred_cells <- colnames(pred_expression)
mse_mat <- matrix(,nrow = length(ref_cells), ncol = length(pred_cells))
rownames(mse_mat) <- ref_cells; colnames(mse_mat) <- pred_cells
for (i in 1:length(ref_cells)) {
  for (j in 1:length(pred_cells)) {
    mse_mat[i,j] <- sum((ref_expression[,i] - pred_expression[,j])^2) / length(gene_list)
  }
}
mse_mat

