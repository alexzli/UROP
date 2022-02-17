datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))

RCTD <- readRDS(file.path(datadir, '/Data/postRCTD_spacexr.rds'))
puck <- RCTD@originalSpatialRNA
RCTD_list <- run.unsupervised(puck, resolution = 0.6, info_type = 'mean', fit_genes = 'de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list, file.path(datadir, 'Objects/MERFISH_15_25.rds'))
