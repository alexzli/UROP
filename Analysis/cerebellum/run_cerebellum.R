datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))
puck <- readRDS(file.path(datadir, '/Data/puckCropped_cerebellum_slideseq.rds'))

RCTD_list <- run.unsupervised(puck, resolution = 0.5)
saveRDS(RCTD_list, file.path(datadir, 'Objects/cerebellum_10_5.rds'))
