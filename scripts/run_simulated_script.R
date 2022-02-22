datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c( .libPaths(), file.path(datadir, 'Rpcg')))
library(spacexr)
source(file.path(datadir, '/R/algorithm.R'))


puck <- readRDS(file.path(datadir, 'data/sim_puck_variableUMI_small.rds'))
RCTD_list <- run.unsupervised(puck, max_cores = 8, resolution = 0.4, SCT = F, silhouette_cutoff = 0, info_type = 'mean', fit_genes = 'de')
saveRDS(RCTD_list, file.path(datadir, 'data/sim_puck_variableUMI_results_4res_SILHOUETTE.rds'))
