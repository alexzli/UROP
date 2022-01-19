library(spacexr)
library(Matrix)
library(doParallel)
library(caret)

source('../iter_optim.R')
source('../analysis.R')

# Load data
reference <- readRDS('../reference_RCTD_vec.rds')
puck <- readRDS('../puckCropped.rds')

# Fit gene expression
myRCTD <- create.RCTD(puck, reference, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
saveRDS(myRCTD,'../myRCTD.rds')
myRCTD <- readRDS('../myRCTD.rds')

myRCTD <- fit.gene.expression(myRCTD)
saveRDS(myRCTD,'../myRCTDE.rds')
myRCTD <- readRDS('../myRCTDE.rds')

# Fit cell types
RCTD2 <- fit.cell.types(myRCTD)
saveRDS(RCTD2,'../RCTD2.rds')
RCTD2 <- readRDS('../RCTD2.rds')

# Comparing ground truth and fitted cell types
cell_types <- c('Astrocytes', 'Bergmann', 'Endothelial', 'Fibroblast', 'Golgi', 'Granule', 'Lugaro', 'Microglia',
                'MLI1', 'MLI2', 'Oligodendrocytes', 'Polydendrocytes', 'Purkinje', 'UBCs')

true_singlet <- myRCTD@results$results_df
true_singlet <- true_singlet[true_singlet$spot_class == 'singlet' & is.element(true_singlet$first_type, cell_types),] # 7106 true singlets
pred_singlet <- RCTD2@results$results_df
pred_singlet <- pred_singlet[pred_singlet$spot_class == 'singlet',] # 8251 predicted singlets
common_barcode <- intersect(row.names(true_singlet), row.names(pred_singlet)) # 6630 common singlets

true_singlet <- true_singlet[common_barcode,]
pred_singlet <- pred_singlet[common_barcode,]

true_types = unlist(list(true_singlet[,'first_type']))
pred_types = unlist(list(pred_singlet[,'first_type']))

conf_mat <- confusionMatrix(pred_types, factor(true_types, levels = levels(pred_types)))
saveRDS(conf_mat,'../conf_mat.rds')

conf_mat <- readRDS('../conf_mat.rds')

# Comparing gene expression of reference to RCTDE
# expression_de : matrix of differential expression genes for cell in cell_types
# expression_ref : matrix of reference genes for cell in cell_types

cell_types <- c('Astrocytes', 'Bergmann', 'Endothelial', 'Fibroblast', 'Golgi', 'Granule', 'Lugaro', 'Microglia',
                'MLI1', 'MLI2', 'Oligodendrocytes', 'Polydendrocytes', 'Purkinje', 'UBCs')

expression_de <- myRCTD@de_results$gene_fits$mean_val
puck = myRCTD@originalSpatialRNA
cti_renorm <- get_norm_ref(puck, myRCTD@cell_type_info$info[[1]], intersect(rownames(expression_de),
	rownames(myRCTD@cell_type_info$info[[1]])), myRCTD@internal_vars$proportions)
expression_ref <- data.matrix(cti_renorm[,cell_types])
expression_ref <- log(expression_ref)
expression_ref <- na.omit(expression_ref)
#plot(expression_de[,'Granule'],expression_ref[,'Granule'])

# Calculating correlation
#colnames(expression_ref) <- c('Granule_ref', 'Purkinje_ref', 'Astrocytes_ref')
#colnames(expression_de) <- c('Granule_de', 'Purkinje_de', 'Astrocytes_de')
expression_tot <- cbind(expression_ref, expression_de)
#expression_tot <- expression_tot[,c('Granule_ref', 'Granule_de', 'Purkinje_ref', 'Purkinje_de', 'Astrocytes_ref', 'Astrocytes_de')]
expression_tot <- as.data.frame(expression_tot)
expression_tot <- expression_tot[expression_tot$Granule_ref != -Inf & expression_tot$Purkinje_ref != -Inf & expression_tot$Astrocytes_ref != -Inf,]
correlation <- cor(expression_tot)
correlation <- as.matrix(correlation)

coul <- pals::brewer.blues(20)[2:20]
heatmap(correlation, symm = TRUE, col = coul, main="Expression Correlation", margins=c(10,10))
