datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c( .libPaths(), file.path(datadir, 'Rpcg')))
library(spacexr)
source(file.path(datadir, '/R/ITERATIVE_OPTIMIZATION.R')

puck <- readRDS(file.path(datadir, '/data/sim_puck.rds'))
RCTD_list <- run.unsupervised(puck, max_cores = 8, info_type = 'de')

saveRDS(RCTDlist, file.path(datadir, '/data/sim_puck_results_de_init.rds'))