library(STdeconvolve)
source('../UROP/R/analysis.R')

# generate STdeconvolve predictions
puck <- readRDS('../UROP/objects/puckCropped_cerebellum_slideseq.rds')
pos <- puck@coords
cd <- puck@counts
counts <- cleanCounts(cd, min.lib.size = 100)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(8, 8, by = 1))
optLDA <- optimalModel(models = ldas, opt = "min")
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
saveRDS(deconProp,'../UROP/objects/deconProp_cerebellum_slideseq.rds')

# function to generate cell assignment table with RCTD and deconvolve results
compare_deconvolve <- function(RCTD, deconProp, verbose = FALSE, singlet_threshold = 0) {
  deconProp <- deconProp[apply(deconProp, 1, max) >= singlet_threshold,]
  deconProp <- as.matrix(apply(deconProp, 1, which.max))
  truth_singlet <- RCTD@results$results_df
  truth_singlet <- truth_singlet[truth_singlet$spot_class == 'singlet',]
  common_barcode <- intersect(row.names(truth_singlet), row.names(deconProp))
  if (verbose) {
    message(paste(dim(truth_singlet)[1], 'ground truth singlets'))
    message(paste(dim(deconProp)[1], 'predicted singlets'))
    message(paste(length(common_barcode), 'common singlets'))
  }
  truth_singlet <- truth_singlet[common_barcode,]
  pred_singlet <- as.matrix(deconProp[common_barcode,])
  colnames(pred_singlet) <- 'first_type'
  truth_types = unlist(list(truth_singlet[,'first_type']))
  pred_types = unlist(list(pred_singlet[,'first_type']))
  table(truth_types, pred_types)
}

# load data
deconProp <- readRDS('../UROP/objects/deconProp_cerebellum_slideseq.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDlist <- readRDS('../UROP/objects/RCTD_list_clustered.rds')
RCTDpred <- RCTDlist[[length(RCTDlist)]]

mytable <- compare_deconvolve(RCTDtruth, deconProp, verbose = TRUE, singlet_threshold = 0.45)
mytable
cell.heatmap(mytable = mytable)
cond.entropy(mytable = mytable)
class.accuracy(mytable = mytable)

cell.assignment.table(RCTDtruth, RCTDpred, verbose = TRUE)
cell.heatmap(RCTDtruth, RCTDpred)
cond.entropy(RCTDtruth, RCTDpred)
class.accuracy(RCTDtruth, RCTDpred)
