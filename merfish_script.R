datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c( .libPaths(), file.path(datadir, 'Rpcg')))
library(spacexr)
library(Matrix)
source(file.path(datadir, '/R/ITERATIVE_OPTIMIZATION.R')

puck <- readRDS(file.path(datadir, '/data/MERFISH_puck_20um2.rds'))
RCTD_list <- run.unsupervised(puck, resolution = 0.18, info_type = 'mean', fit_genes = 'mean', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,file.path(datadir, '/data/MERFISH20um2_results_mean_fit_fullgenes.rds'))