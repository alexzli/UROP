library(STdeconvolve)
source('../UROP/R/analysis.R')

# generate STdeconvolve predictions
puck <- readRDS('../UROP/data/puckCropped_cerebellum_slideseq.rds')
pos <- puck@coords
cd <- puck@counts
counts <- cleanCounts(cd, min.lib.size = 100)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = 12)
optLDA <- optimalModel(models = ldas, opt = "min")
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
saveRDS(results,'../UROP/objects/stdeconvolve_cerebellum_slideseq_10res.rds')

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
  return(table(truth_types, pred_types))
}

# load data
results <- readRDS('../UROP/objects/stdeconvolve_cerebellum_slideseq_10res.rds')
deconProp <- results$theta
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDlist <- readRDS('../UROP/results/cerebellum_slideseq_results_mean.rds')
RCTDpred <- RCTDlist[[length(RCTDlist)]]

mytable <- compare_deconvolve(RCTDtruth, deconProp, verbose = TRUE, singlet_threshold = 0.35)
mytable
cell.heatmap(mytable = mytable)
cond.entropy(mytable = mytable)
class.accuracy(mytable = mytable)

mytable <- cell.table(RCTDtruth, RCTDpred, verbose = TRUE)
cell.heatmap(RCTDtruth, RCTDpred)
cond.entropy(RCTDtruth, RCTDpred)
class.accuracy(RCTDtruth, RCTDpred)


mytable <- rbind(mytable, mytable['MLI1',] + mytable['MLI2',])
rownames(mytable)[length(rownames(mytable))] <- 'MLI'
mytable <- mytable[-c(which(rownames(mytable) == 'MLI1'), which(rownames(mytable) == 'MLI2')),]
assignment.accuracy(mytable)
assignment.accuracy.plot(mytable)


##compare gene expression of stdeconvolve with ground truth with iterative alg
# preprocessing
results <- readRDS('../UROP/objects/stdeconvolve_cerebellum_slideseq_13res.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDlist <- readRDS('../UROP/objects/RCTD_list_clustered_13res.rds')
ref_info <- RCTDtruth@cell_type_info$renorm[[1]]
pred_info <- t(results$beta)/1000
gene_list <- intersect(rownames(ref_info), rownames(pred_info))
length(gene_list)
ref_info <- data.matrix(ref_info[gene_list,])
pred_info <- data.matrix(pred_info[gene_list,])

# stdeconvolve plots
gene_expression <- list(reference = ref_info, prediction = pred_info)
mytable <- compare_deconvolve(RCTDtruth, results$theta, verbose = TRUE, singlet_threshold = 0.3)
#gene.aggregate.mse(gene_expression, mytable)
rank <- aggregate.gene.rank(gene_expression, mytable)
gene.bin.plot(rank, rank_lim=100, bins=20)
gene.density.plot(rank, rank_lim=100)

# iterative algorithm plots
gene_expression <- gene.expression(RCTDtruth, RCTDlist[[length(RCTDlist)]], gene_list = gene_list)
mytable <- cell.table(RCTDtruth, RCTDlist[[length(RCTDlist)]], verbose = TRUE)
#gene.aggregate.mse(gene_expression, mytable)
rank <- aggregate.gene.rank(gene_expression, mytable)
gene.bin.plot(rank, rank_lim=100, bins=20)
gene.density.plot(rank, rank_lim=100)

