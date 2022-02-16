source('../UROP/R/ITERATIVE_OPTIMIZATION.R')

# Old code for running the iterative algorithm on various datasets. More current
# code running has been moved to the datasets' respective folders.


# new MERFISH puck
RCTD <- readRDS('../UROP/data/postRCTD_spacexr.rds')
puck <- RCTD@originalSpatialRNA
RCTD_list <- run.unsupervised(puck, resolution = 0.6, info_type = 'mean', fit_genes = 'de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,'../UROP/results/NEW_MERFISH.rds')

# simulated puck
puck <- readRDS('../UROP/objects/sim_puck.rds')
RCTD_list <- run.unsupervised(puck, resolution = 0.15, info_type = 'mean', fit_genes = 'singlet de')
saveRDS(RCTD_list,'../UROP/results/sim_puck_results_singletde_fit.rds')

# simulated MERFISH puck
puck <- readRDS('../UROP/objects/MERFISH_puck_10um2.rds')
RCTD_list <- run.unsupervised(puck, resolution = 0.15, info_type = 'mean', fit_genes = 'de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,'../UROP/results/MERFISH10um2_results_de_fit_NEWCUTOFFS.rds')

puck <- readRDS('../UROP/objects/MERFISH_puck_20um2.rds')
RCTD_list <- run.unsupervised(puck, resolution = 0.08, SCT = T, info_type = 'mean', fit_genes = 'de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,'../UROP/results/MERFISH20um2_results_de_fit_2_10_cutoff.rds')

puck <- readRDS('../UROP/objects/MERFISH_puck_10um2.rds')
RCTD_list <- run.unsupervised(puck, resolution = 0.15, info_type = 'mean', fit_genes = 'mean', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,'../UROP/results/MERFISH10um2_results_mean_fit_NEWCUTOFFS.rds')

puck <- readRDS('../UROP/objects/MERFISH_puck_20um2.rds')
RCTD_list <- run.unsupervised(puck, resolution = 0.18, info_type = 'mean', fit_genes = 'de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,'../UROP/results/MERFISH20um2_results_mean_fit_2_10.rds')

puck <- readRDS('../UROP/objects/MERFISH_puck_10um2.rds')
RCTD_list <- run.unsupervised(puck, resolution = 0.15, info_type = 'mean', fit_genes = 'de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,'../UROP/results/MERFISH10um2_results_singletde_fit_fullgenes.rds')

puck <- readRDS('../UROP/objects/MERFISH_puck_20um2.rds')
RCTD_list <- run.unsupervised(puck, resolution = 0.18, info_type = 'mean', fit_genes = 'singlet de', gene_list = puck@counts@Dimnames[[1]])
saveRDS(RCTD_list,'../UROP/results/MERFISH20um2_results_singletde_fit_fullgenes.rds')

# cerebellum slide-seq
puck <- readRDS('../UROP/data/puckCropped_cerebellum_slideseq.rds')

RCTD_list <- run.unsupervised(puck, info_type = 'mean')
saveRDS(RCTD_list,'../UROP/results/cerebellum_slideseq_results_mean.rds')

RCTD_list <- run.unsupervised(puck, resolution = 1, info_type = 'mean', fit_genes = 'de')
saveRDS(RCTD_list,'../UROP/results/cerebellum_slideseq_results_de.rds')

RCTD_truth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTD_list <- run.semisupervised(RCTD_truth, keep_gene_list = T)
saveRDS(RCTD_list,'../UROP/results/cerebellum_slideseq_truth_list.rds')

# cancer dataset
myRCTD <- readRDS('../UROP/data/myRCTD_201014_03.rds')
results <- read.csv(file = '../UROP/data/results_processed_slideseq_data_2020-12-12_Puck_201014_03_2020-12-12_Puck_201014_03_cell_types.csv')
cancer_results <- results[results$HRS_location == TRUE,]
results <- myRCTD@results$results_df
results$first_type = factor(results$first_type, levels=append(levels(results$first_type), 'Cancer_cells'))
results$second_type = factor(results$second_type, levels=append(levels(results$second_type), 'Cancer_cells'))
for (barcode in cancer_results$barcodes) {
  results[barcode, 'spot_class'] = 'singlet'
  results[barcode,'first_type'] = 'Cancer_cells'
}
myRCTD@results$results_df <- results
myRCTD@cell_type_info$info[[2]] <- append(myRCTD@cell_type_info$info[[2]], 'Cancer_cells')
myRCTD@cell_type_info$renorm[[2]] <- append(myRCTD@cell_type_info$renorm[[2]], 'Cancer_cells')
myRCTD@cell_type_info$info[[3]] <- myRCTD@cell_type_info$info[[3]] + 1
myRCTD@cell_type_info$renorm[[3]] <- myRCTD@cell_type_info$renorm[[3]] + 1
class_df <- myRCTD@internal_vars$class_df
class_df['Cancer_cells', 'class'] = 'Cancer_cells'
myRCTD@internal_vars$class_df <- class_df
myRCTD@config$max_cores <- 4

RCTD_list <- run.semisupervised(myRCTD, keep_gene_list = T)
saveRDS(RCTD_list,'../UROP/results/cancer_results_keepgenes.rds')

#cell_type_info <- fit.gene.expression(myRCTD)
#RCTD <- create.RCTD.noref(myRCTD@originalSpatialRNA, max_cores = 4, cell_type_info = cell_type_info, 
#                            gene_list_reg = myRCTD@internal_vars$gene_list_reg, class_df = myRCTD@internal_vars$class_df)
#iter.optim(RCTD)

RCTD_list <- run.semisupervised(myRCTD, keep_gene_list = F)
saveRDS(RCTD_list,'../UROP/results/cancer_results_diffgenes.rds')
