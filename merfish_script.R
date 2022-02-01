datadir <- '/broad/thechenlab/Alex/UROP/'
.libPaths( c( .libPaths(), file.path(datadir, 'Rpcg')))
library(spacexr)
library(Matrix)
source(file.path(datadir, '/R/optimization.R'))
source(file.path(datadir, '/R/initialization.R'))
source(file.path(datadir, '/R/clustering.R'))
source(file.path(datadir, '/R/analysis.R'))


simFN7 <- readRDS(file.path(datadir, '/data/simFN7_10um2.rds'))
puck <- SpatialRNA(simFN7$cellCounts[,1:2], t(simFN7$sim), require_int = FALSE)
assignments <- gen.clusters(puck, resolution = 0.3, SCT = F)
cell_type_info <- cell_type_info_from_assignments(puck, assignments)
myRCTD <- create.RCTD.noref(puck, cell_type_info = cell_type_info)
myRCTD@config$max_cores <- 32

RCTDlist <- run.iter.optim(myRCTD, cell_types_assigned = FALSE, max_n_iter = 50, CELL_MIN_INSTANCE = 0, constant_genes = TRUE)
saveRDS(RCTDlist, file.path(datadir, '/data/MERFISH_10um2_results.rds'))