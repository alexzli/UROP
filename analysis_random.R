source('../UROP/R/optimization.R')
source('../UROP/R/initialization.R')
source('../UROP/R/clustering.R')
source('../UROP/R/analysis.R')

# generate results
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
gene_list_reg <- RCTDtruth@internal_vars$gene_list_reg
gene_list_bulk <- RCTDtruth@internal_vars$gene_list_bulk

puck <- readRDS('../UROP/objects/puckCropped_cerebellum_slideseq.rds')
cell_type_info <- gen.random.info(puck, n_cell_types = 8)
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info, gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk)
RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, max_n_iter = 20, CELL_MIN_INSTANCE = 0, constant_genes = TRUE)
saveRDS(RCTDlist,'../UROP/objects/RCTD_list_random.rds')

# load results
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDlist <- readRDS('../UROP/objects/RCTD_list_random.rds')


# Cell type assignments compare to initial (testing deviation from initialization point)
for (i in 2:length(RCTDlist)) {
  print(cell.confusion.mat(RCTDlist[[1]], RCTDlist[[i]])$overall['Accuracy'])
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
