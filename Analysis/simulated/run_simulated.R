datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))

puck <- readRDS(file.path(datadir, 'Data/sim_puck.rds'))
RCTD_list <- run.unsupervised(puck, resolution = 0.15, info_type = 'mean', fit_genes = 'de')
saveRDS(RCTD_list, file.path(datadir, 'Objects/sim_puck_results_de_fit.rds'))