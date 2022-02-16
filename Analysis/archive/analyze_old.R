datadir <- '/UROP/'
source(file.path(datadir, '/R/analysis.R'))

# Old code for analyzing results from various datasets. More current
# code running has been moved to the datasets' respective folders.


# generate truth types
reference <- readRDS(file.path(datadir, '/data/reference_RCTD_vec.rds'))
RCTDtruth <- create.RCTD(puck, reference, max_cores = 4)
RCTDtruth <- run.RCTD(RCTDtruth, doublet_mode = 'doublet')
saveRDS(RCTDtruth,file.path(datadir, '/objects/RCTD_cerebellum_slideseq.rds'))

# generate results (cell type info is input instead of cell types)
puck <- readRDS('../UROP/data/puckCropped_cerebellum_slideseq.rds')
assignments <- gen.clusters(puck, resolution = 1.3, SCT = F)
cell_type_info <- cell_type_info_from_assignments(puck, assignments)
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, max_n_iter = 20, CELL_MIN_INSTANCE = 25, constant_genes = TRUE)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_clustered_13res_genethresh0.rds')

# generate results (cell types directly assigned)
puck <- readRDS('../UROP/objects/puckCropped_cerebellum_slideseq.rds')
assignments <- gen.clusters(puck, resolution = 1.3, SCT = F)
RCTDcluster <- create.RCTD.noref(puck)
RCTDcluster <- set_internal_vars(RCTDcluster)
RCTDcluster <- assign.cell.types(RCTDcluster, assignments)
RCTDcluster <- fit.gene.expression(RCTDcluster, levels(RCTDcluster@results$results_df$first_type))
cell_type_info <- list(as.data.frame(exp(RCTDcluster@de_results$gene_fits$mean_val)), levels(RCTDcluster@results$results_df$first_type), length(levels(RCTDcluster@results$results_df$first_type)))
cell_type_info <- list(info = cell_type_info, renorm = cell_type_info)
saveRDS(RCTDcluster,'../UROP/objects/RCTD_cluster_initial_assignments_23res.rds')
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, max_n_iter = 10, CELL_MIN_INSTANCE = 25, constant_genes = TRUE)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_de_generated_info.rds')


# load results
RCTDlist <- readRDS('../UROP/objects/RCTD_list_clustered_13res.rds')
RCTDpred <- RCTDlist[[length(RCTDlist)]]
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDcluster <- readRDS('../UROP/objects/RCTD_cluster_initial_assignments_3res.rds')


# Cell type assignments compare to initial (testing deviation from initialization point)
for (i in 1:length(RCTDlist)) {
  print(cell.confusion.mat(RCTDcluster, RCTDlist[[i]])$overall['Accuracy'])
}

# Cell type assignments compare consecutive (testing convergence)
for (i in 1:(length(RCTD_list)-1)) {
  print(cell.confusion.mat(RCTD_list[[i]], RCTD_list[[i+1]])$overall['Accuracy'])
}

# Cell type assignments compare to ground truth (conditional entropy)
for (i in 1:length(RCTDlist)) {
  print(cond.entropy(RCTDtruth, RCTDlist[[i]]))
}

# Cell type assignments compare to ground truth (classification accuracy)
for (i in 1:length(RCTDlist)) {
  print(class.accuracy(RCTDtruth, RCTDlist[[i]]), verbose = TRUE)
}

# Cell type assignments compare to ground truth (table)
for (i in 1:length(RCTDlist)) {
  print(cell.assignment.table(RCTDtruth, RCTDlist[[i]], verbose = TRUE))
}

# Gene expression MSE compare to initial (testing deviation from initialization point)
for (i in 2:length(RCTD_list)) {
  print(gene.mse(RCTD_list[[1]], RCTD_list[[i]]))
}

# Gene expression MSE compare consecutive (testing convergence)
for (i in 1:(length(RCTDlist)-1)) {
  print(gene.mse(RCTDlist[[i]], RCTDlist[[i+1]]))
}

# Pixel assignment accuracy compare consecutive (testing convergence)
for (i in 1:(length(RCTD_list)-1)) {
  print(assignment_accuracy(RCTD_list[[i]], RCTD_list[[i+1]]))
}

# spot type distribution over iterations
for (i in 1:length(RCTD_list)) {
  print(table(RCTD_list[[i]]@results$results_df$spot_class))
}

# gene correlation over iterations
for (i in 1:length(RCTDlist)) {
  print(gene.heatmap(RCTDtruth, RCTDlist[[i]]))
}

# compare gene expression ranks
gene_expression <- gene.expression(RCTDtruth, RCTDlist[[4]])
mytable <- cell.table(RCTDtruth, RCTDlist[[4]])
rank <- aggregate.gene.rank(gene_expression, mytable)
gene.aggregate.mse(gene_expression, mytable)
gene.bin.plot(rank, rank_lim=1000, bins=20)
gene.density.plot(rank, rank_lim=1000)

for (i in 1:length(RCTDlist)) {
  gene_expression <- gene.expression(RCTDtruth, RCTDlist[[i]])
  mytable <- cell.table(RCTDtruth, RCTDlist[[i]])
  print(gene.aggregate.mse(gene_expression, mytable))
}

gene_expression <- gene.expression(RCTDtruth, RCTDlist[[length(RCTDlist)]])
mytable <- cell.table(RCTDtruth, RCTDlist[[length(RCTDlist)]])
cluster_assignments <- apply(mytable, 2, function(x) which(x == max(x)))
rank <- gene.expression.rank(gene_expression, 16, 4)

# load results
RCTDlist <- readRDS('../UROP/objects/RCTD_list_clustered_13res.rds')
RCTDpred <- RCTDlist[[length(RCTDlist)]]
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDcluster <- readRDS('../UROP/objects/RCTD_cluster_initial_assignments_13res.rds')

# COMBINE MLI1 AND MLI2
mytable <- cell.table(RCTDtruth, RCTDcluster)
mytable <- rbind(mytable, mytable['MLI1',] + mytable['MLI2',])
rownames(mytable)[length(rownames(mytable))] <- 'MLI'
mytable <- mytable[-c(which(rownames(mytable) == 'MLI1'), which(rownames(mytable) == 'MLI2')),]
assignment.accuracy(mytable)
assignment.accuracy.plot(mytable)

for (i in 1:length(RCTDlist)) {
  mytable <- cell.table(RCTDtruth, RCTDlist[[i]])
  mytable <- rbind(mytable, mytable['MLI1',] + mytable['MLI2',])
  rownames(mytable)[length(rownames(mytable))] <- 'MLI'
  mytable <- mytable[-c(which(rownames(mytable) == 'MLI1'), which(rownames(mytable) == 'MLI2')),]
  print(assignment.accuracy(mytable))
}
