source('../UROP/R/optimization.R')
source('../UROP/R/initialization.R')
source('../UROP/R/clustering.R')
source('../UROP/R/analysis.R')

# generate results (cell type info is input instead of cell types)
puck <- readRDS('../UROP/objects/puckCropped_cerebellum_slideseq.rds')
assignments <- gen.clusters(puck, SCT = F)
cell_type_info <- cell_type_info_from_assignments(puck, assignments)
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, max_n_iter = 10, CELL_MIN_INSTANCE = 25, constant_genes = TRUE)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_clustered.rds')

reference <- readRDS('../UROP/objects/reference_RCTD_vec.rds')
RCTDtruth <- create.RCTD(puck, reference, max_cores = 4)
RCTDtruth <- run.RCTD(RCTDtruth, doublet_mode = 'doublet')
saveRDS(RCTDtruth,'../UROP/objects/RCTD_cerebellum_slideseq.rds')

# load results
RCTDlist <- readRDS('../UROP/objects/RCTD_list_clustered.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')


# generate results (cell types directly assigned)
puck <- readRDS('../UROP/objects/puckCropped_cerebellum_slideseq.rds')
assignments <- gen.clusters(puck, SCT = F)
myRCTD <- create.RCTD.noref(puck)
myRCTD <- set_internal_vars(myRCTD)
myRCTD <- assign.cell.types(myRCTD, assignments)
myRCTD <- run.iter.optim(myRCTD, cell_types_assigned = TRUE, max_n_iter = 10, CELL_MIN_INSTANCE = 25, constant_genes = TRUE)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_de_generated_info.rds')
# load results
RCTDlist <- readRDS('../UROP/objects/RCTD_list_de_generated_info.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')


# Cell type assignments compare to initial (testing deviation from initialization point)
for (i in 2:length(RCTDlist)) {
  print(cell.confusion.mat(RCTDlist[[1]], RCTDlist[[i]])$overall['Accuracy'])
}

# Cell type assignments compare consecutive (testing convergence)
for (i in 1:(length(RCTDlist)-1) {
  print(cell.confusion.mat(RCTDlist[[i]], RCTDlist[[i+1]])$overall['Accuracy'])
}

# Cell type assignments compare to ground truth (conditional entropy)
for (i in 1:length(RCTDlist)) {
  print(cond.entropy(RCTDtruth, RCTDlist[[i]]))
}

# Cell type assignments compare to ground truth (classification accuracy)
for (i in 1:length(RCTDlist)) {
  print(class.accuracy(RCTDtruth, RCTDlist[[i]], verbose = TRUE))
}

# Cell type assignments compare to ground truth (table)
for (i in 1:length(RCTDlist)) {
  print(cell.assignment.table(RCTDtruth, RCTDlist[[i]]))
}

# Gene expression MSE compare to initial (testing deviation from initialization point)
for (i in 2:length(RCTDlist)) {
  print(gene.mse(RCTDlist[[1]], RCTDlist[[i]]))
}

# Gene expression MSE compare consecutive (testing convergence)
for (i in 1:(length(RCTDlist)-1)) {
  print(gene.mse(RCTDlist[[i]], RCTDlist[[i+1]]))
}
