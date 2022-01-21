source('../UROP/R/optimization.R')
source('../UROP/R/initialization.R')
source('../UROP/R/clustering.R')
source('../UROP/R/analysis.R')

# generate results (cell type info is input instead of cell types)
puck <- readRDS('../UROP/objects/puckCropped_cerebellum_slideseq.rds')
assignments <- gen.clusters(puck, resolution = 2.3, SCT = F)
cell_type_info <- cell_type_info_from_assignments(puck, assignments)
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, max_n_iter = 20, CELL_MIN_INSTANCE = 25, constant_genes = TRUE)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_clustered_23res.rds')

reference <- readRDS('../UROP/objects/reference_RCTD_vec.rds')
RCTDtruth <- create.RCTD(puck, reference, max_cores = 4)
RCTDtruth <- run.RCTD(RCTDtruth, doublet_mode = 'doublet')
saveRDS(RCTDtruth,'../UROP/objects/RCTD_cerebellum_slideseq.rds')

# load results
RCTDlist <- readRDS('../UROP/objects/RCTD_list_clustered_3res.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDcluster <- readRDS('../UROP/objects/RCTD_cluster_initial_assignments_3res.rds')


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
RCTDlist2 <- readRDS('../UROP/objects/RCTD_list_de_generated_info.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')


# Cell type assignments compare to initial (testing deviation from initialization point)
for (i in 1:length(RCTDlist)) {
  print(cell.confusion.mat(RCTDcluster, RCTDlist[[i]])$overall['Accuracy'])
}

# Cell type assignments compare consecutive (testing convergence)
for (i in 1:(length(RCTDlist)-1)) {
  print(cell.confusion.mat(RCTDlist[[i]], RCTDlist[[i+1]])$overall['Accuracy'])
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
for (i in 2:length(RCTDlist)) {
  print(gene.mse(RCTDlist[[1]], RCTDlist[[i]]))
}

# Gene expression MSE compare consecutive (testing convergence)
for (i in 1:(length(RCTDlist)-1)) {
  print(gene.mse(RCTDlist[[i]], RCTDlist[[i+1]]))
}

