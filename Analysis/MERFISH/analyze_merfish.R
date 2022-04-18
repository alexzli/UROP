datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))
source(file.path(datadir, 'R/analysis.R'))

RCTD_list <- readRDS(file.path(datadir, 'Objects/MERFISH_15_25_SCT.rds'))
RCTDpred <- RCTD_list[[length(RCTD_list)]]

## INITIAL DATA ANALYSIS

# Spatially plot predicted cells
plot.cell.types(RCTDpred)
# Gene expression MSE compare consecutive (testing convergence)
for (i in 1:(length(RCTD_list)-1)) {
  print(gene.mse(RCTD_list[[i]], RCTD_list[[i+1]]))
}
# Pixel assignment accuracy compare consecutive (testing convergence)
for (i in 1:(length(RCTD_list)-1)) {
  print(assignment_accuracy(RCTD_list[[i]], RCTD_list[[i+1]]))
}
# spot type distribution over iterations
for (i in 1:length(RCTD_list)) {
  print(table(RCTD_list[[i]]@results$results_df$spot_class))
}


## COMPARE TO GROUND TRUTH

# load truth data
myRCTD <- readRDS(file.path(datadir, '/Data/postRCTD_spacexr.rds'))
moffit <- read.csv2(file = file.path(datadir, "/Data/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv"),
                    sep = ",")
ids <- colnames(myRCTD@originalSpatialRNA@counts)
truth <- moffit[moffit$Cell_ID %in% ids,][,c(1, 8)]
rownames(truth) <- truth$Cell_ID
truth$Cell_ID <- NULL

# singlet analysis
library(tm)
singlet_table <- function(RCTD) {
  singlet_ids <- rownames(RCTD@results$results_df[RCTD@results$results_df$spot_class == 'singlet',])
  pred_types <- RCTD@results$results_df[singlet_ids,'first_type']
  true_types <- truth[singlet_ids,'Cell_class']
  true_types <- removeNumbers(true_types)
  mytable <- table(true_types, pred_types) 
  return(mytable)
}
mytable <- singlet_table(RCTDpred)
assignment.accuracy(mytable)
assignment.accuracy.plot(mytable)
mytable

for (i in 1:length(RCTD_list)) {
  mytable <- singlet_table(RCTD_list[[i]])
  print(assignment.accuracy(mytable))
}


# excitatory and inhibitory neurons across iterations
for (i in 1:length(RCTD_list)) {
  mytable <- singlet_table(RCTD_list[[i]])
  print(paste('iteration', i, ':', assignment.accuracy(mytable)))
  print(mytable[c('Excitatory', 'Inhibitory'),])
}

# analyze singlet scores
high_scores <- RCTDpred@results$results_df[RCTDpred@results$results_df$singlet_score - RCTDpred@results$results_df$min_score > 40,]
truth[rownames(high_scores),] # ground truth cell types of spots with very high singlet - min score


# compare with clustering results
myRCTD <- readRDS(file.path(datadir, '/Data/postRCTD_spacexr.rds'))
assignments <- gen.clusters(myRCTD@originalSpatialRNA, resolution = 0.15, silhouette_cutoff = 0)
RCTDcluster <- create.RCTD.noref(myRCTD@originalSpatialRNA)
RCTDpred2 <- assign.cell.types(RCTDcluster, assignments)

mytable2 <- singlet_table(RCTDpred2)
assignment.accuracy(mytable2)
assignment.accuracy.plot(mytable2)
mytable2

