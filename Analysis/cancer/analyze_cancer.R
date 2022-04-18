datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))
source(file.path(datadir, 'R/analysis.R'))

RCTDpred <- readRDS(file.path(datadir, 'Objects/cancer_final.rds'))
plot_cell_types(RCTDpred)

## IMMUNE CELL SUBTYPES

# CD4 T cells subtypes
CD4_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet' &
                                              RCTDpred@results$results_df$first_type == 'CD4_T_cells',])
CD4_puck <- restrict_puck(RCTDpred@originalSpatialRNA, CD4_barcodes)
CD4_RCTD <- run.unsupervised(
  CD4_puck,
  resolution = 0.65,
  doublet_mode = 'full',
  max_iter = 200
)
saveRDS(CD4_RCTD, file.path(datadir, 'Objects/cancer_CD4subtypes.rds'))
CD4_RCTD <- readRDS(file.path(datadir, 'Objects/cancer_CD4subtypes.rds'))
plot_puck_continuous(CD4_puck, CD4_barcodes, CD4_RCTD@results$weights[,4], my_pal = pals::brewer.blues(20)[2:20])

# CD8T cells subtypes
CD8_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet' &
                                                       RCTDpred@results$results_df$first_type == 'CD8_T_cells',])
CD8_puck <- restrict_puck(RCTDpred@originalSpatialRNA, CD8_barcodes)
CD8_RCTD <- run.unsupervised(
  CD8_puck,
  resolution = 0.75,
  doublet_mode = 'full',
  max_iter = 200
)
saveRDS(CD8_RCTD, file.path(datadir, 'Objects/cancer_CD8subtypes.rds'))
CD8_RCTD <- readRDS(file.path(datadir, 'Objects/cancer_CD8subtypes.rds'))
plot_puck_continuous(CD8_puck, CD8_barcodes, CD8_RCTD@results$weights[,3], my_pal = pals::brewer.blues(20)[2:20])

# Monocytes/macrophaages cells subtypes
MM_barcodes <- rownames(RCTDpred@results$results_df[RCTDpred@results$results_df$spot_class == 'singlet' &
                                                       RCTDpred@results$results_df$first_type == 'Monocytes_macrophages',])
MM_puck <- restrict_puck(RCTDpred@originalSpatialRNA, MM_barcodes)
MM_RCTD <- run.unsupervised(
  MM_puck,
  resolution = 0.7,
  doublet_mode = 'full',
  max_iter = 200
)
saveRDS(MM_RCTD, file.path(datadir, 'Objects/cancer_MMsubtypes.rds'))
MM_RCTD <- readRDS(file.path(datadir, 'Objects/cancer_MMsubtypes.rds'))
plot_puck_continuous(MM_puck, MM_barcodes, MM_RCTD@results$weights[,1], my_pal = pals::brewer.blues(20)[2:20])
