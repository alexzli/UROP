datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))

puck <- readRDS(file.path(datadir, '/Objects/MERFISH_puck_20um2.rds')) # 9 ground truth types, 0.15 res for 20 um2
gene_list <- puck@counts@Dimnames[[1]]
RCTD_list <- run.unsupervised(
  puck,
  resolution = 0.15,
  gene_list = gene_list
)
saveRDS(RCTD_list, file.path(datadir, 'Objects/MERFISH_20um2_10_5.rds'))

