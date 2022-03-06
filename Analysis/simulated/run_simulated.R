datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))

puck <- readRDS(file.path(datadir, 'objects/sim_puck_variableUMI_small.rds'))
RCTD_list <- run.unsupervised(puck, resolution = 0.3, SCT = F, silhouette_cutoff = 0, info_type = 'mean', fit_genes = 'de')
saveRDS(RCTD_list, file.path(datadir, 'objects/sim_puck_variableUMI_results_3res_SILHOUETTE.rds'))
