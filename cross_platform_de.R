library(spacexr)
library(Matrix)
library(doParallel)

source('../iter_optim.R')
source('../analysis.R')

# Load data
train_ref <- readRDS('../reference_RCTD_vec.rds') # reference data
test_ref <- readRDS('../cerebellum_reference_dropviz_spacexr.rds') # data to learn
counts <- test_ref@counts
nUMI <- test_ref@nUMI
puck <- SpatialRNA(NULL, counts, nUMI, use_fake_coords = TRUE)

# Fit gene expression
myRCTD <- create.RCTD(spatialRNA, reference, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
saveRDS(myRCTD,'../cross_platform_RCTD.rds')
myRCTD <- readRDS('../cross_platform_RCTD.rds')

myRCTD <- fit.gene.expression(myRCTD)
saveRDS(myRCTD,'../cross_platform_RCTDE.rds')
myRCTD <- readRDS('../cross_platform_RCTDE.rds')

# Comparing gene expression of testing to training datasets
# test_expression : matrix of cell expression profiles for testing dataset (cerebellum scRNA-seq) after RCTDE
# train_expression : matrix of cell expression profiles for training dataset (cerebellum snRNA-seq)
cell_types <- c('Astrocytes', 'Bergmann', 'Endothelial', 'Fibroblast', 'Granule',
	'MLI1', 'Oligodendrocytes')
test_expression <- myRCTD@de_results$gene_fits$mean_val
puck = myRCTD@originalSpatialRNA
cti_renorm <- get_norm_ref(puck, myRCTD@cell_type_info$info[[1]], intersect(rownames(test_expression),
	rownames(myRCTD@cell_type_info$info[[1]])), myRCTD@internal_vars$proportions)
train_expression <- data.matrix(cti_renorm[,cell_types])
train_expression <- log(train_expression)
train_expression <- na.omit(train_expression)
test_expression <- test_expression[rownames(train_expression),]
#plot(test_expression[,'Granule'],train_expression[,'Granule'])

# Calculating correlation
colnames(train_expression) <- c('Astrocytes_train', 'Bergmann_train', 'Endothelial_train', 'Fibroblast_train', 'Granule_train',
	'MLI1_train', 'Oligodendrocytes_train')
colnames(test_expression) <- c('Astrocytes_test', 'Bergmann_test', 'Endothelial_test', 'Fibroblast_test', 'Granule_test',
	'MLI1_test', 'Oligodendrocytes_test')
expression_tot <- cbind(train_expression, test_expression)
#expression_tot <- expression_tot[,c('Astrocytes_train', 'Astrocytes_test', 'Bergmann_train', 'Bergmann_test', 'Endothelial_train', 'Endothelial_test',
# 'Fibroblast_train', 'Fibroblast_test', 'Granule_train', 'Granule_test', 'MLI1_train', 'MLI1_test', 'Oligodendrocytes_train', 'Oligodendrocytes_test')]
expression_tot <- as.data.frame(expression_tot)
expression_tot <- expression_tot[expression_tot$Astrocytes_train != -Inf & expression_tot$Bergmann_train != -Inf & expression_tot$Endothelial_train != -Inf & expression_tot$Fibroblast_train != -Inf & expression_tot$Granule_train != -Inf & expression_tot$MLI1_train != -Inf & expression_tot$Oligodendrocytes_train != -Inf,]
correlation <- cor(expression_tot)
correlation <- as.matrix(correlation)

test_correlation <- correlation[8:14,8:14]
train_correlation <- correlation[1:7,1:7]
cross_correlation <- correlation[1:7,8:14]

coul <- pals::brewer.blues(20)[2:20]
heatmap(cross_correlation, col = coul, symm = TRUE, main="Transformed Testing vs. Training Dataset Correlation", margins=c(10,10))


# Comparing gene expression of original datasets
# test_original : matrix of cell expression profiles for testing dataset (cerebellum scRNA-seq) before RCTDE
test_ref <- readRDS('../cerebellum_reference_dropviz_spacexr.rds') # data to learn
counts <- test_ref@counts
nUMI <- test_ref@nUMI
puck <- SpatialRNA(NULL, counts, nUMI, use_fake_coords = TRUE)

refRCTD <- create.RCTD(puck, test_ref, max_cores = 4, CELL_MIN_INSTANCE = 0)

cell_types <- c('Astrocytes', 'Bergmann', 'Endothelial', 'Fibroblast', 'Granule',
  'MLI1', 'Oligodendrocytes')
puck = refRCTD@originalSpatialRNA
cti_renorm <- get_norm_ref(puck, refRCTD@cell_type_info$info[[1]], intersect(rownames(test_expression),
  rownames(refRCTD@cell_type_info$info[[1]])), myRCTD@internal_vars$proportions)
test_original <- data.matrix(cti_renorm[,cell_types])
test_original <- log(test_original)
test_original <- na.omit(test_original)

colnames(test_original) <- c('Astrocytes_original', 'Bergmann_original', 'Endothelial_original', 'Fibroblast_original',
  'Granule_original', 'MLI1_original', 'Oligodendrocytes_original')
expression <- cbind(train_expression, test_expression, test_original)
expression <- as.data.frame(expression)
expression <- expression[expression$Astrocytes_train != -Inf & expression$Bergmann_train != -Inf & expression$Endothelial_train != -Inf & expression$Fibroblast_train != -Inf & expression$Granule_train != -Inf & expression$MLI1_train != -Inf & expression$Oligodendrocytes_train != -Inf,]
expression <- expression[expression$Astrocytes_original != -Inf & expression$Bergmann_original != -Inf & expression$Endothelial_original != -Inf & expression$Fibroblast_original != -Inf & expression$Granule_original != -Inf & expression$MLI1_original != -Inf & expression$Oligodendrocytes_original != -Inf,]
correlation <- cor(expression)
correlation <- as.matrix(correlation)

original_test_correlation <- correlation[15:21,15:21]
transformed_test_correlation <- correlation[1:7,15:21]
original_cross_correlation <- correlation[8:14,15:21]

coul <- pals::brewer.blues(20)[2:20]
heatmap(original_cross_correlation, col = coul, symm = TRUE, main="Original Testing vs. Training Dataset Correlation", margins=c(10,10))
