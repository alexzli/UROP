datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c( .libPaths(), file.path(datadir, 'Rpcg')))
library(spacexr)
library(Matrix)
source('/broad/thechenlab/Alex/spacexr/AnalysisCSIDE/helper_functions/de_simulation_helper.R')

reference <- readRDS(file.path(datadir, 'data/reference_RCTD_vec.rds'))
cell_types <- c('Astrocytes', 'Bergmann', 'Fibroblast', 'Golgi', 'Granule', 'MLI1', 'MLI2', 'Oligodendrocytes', 'Polydendrocytes', 'Purkinje')

puck <- generate_sim_puck(cell_types, rownames(reference@counts), reference, trials=30)
saveRDS(puck, file.path(datadir, 'data/sim_puck.rds'))
