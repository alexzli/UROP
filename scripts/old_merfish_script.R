datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c(file.path(datadir, 'Rpcg'), .libPaths()))
source(file.path(datadir, '/R/ITERATIVE_OPTIMIZATION.R'))

puck <- readRDS(file.path(datadir, '/data/MERFISH_puck_20um2.rds'))
RCTD_list <- run.unsupervised(puck, max_cores = 8, resolution = 0.08, SCT = T, info_type = 'mean', fit_genes = 'de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,file.path(datadir, '/data/MERFISH20um2_results_de_fit_2_10_SCT.rds'))