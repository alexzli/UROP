library(spacexr)
library(Matrix)
library(doParallel)

source('../iter_optim_TEST.R')
source('../analysis.R')

# Run on Puck and snRNA-seq reference
myRCTD <- readRDS('../myRCTD.rds')
barcodes <- intersect(names(myRCTD@spatialRNA@nUMI), colnames(myRCTD@spatialRNA@counts))
cell_type_count <- aggregate_cell_types(myRCTD, barcodes, doublet_mode = TRUE)
cell_types <- names(which(cell_type_count >= 25))
RCTD <- fit.gene.expression(myRCTD, cell_types, CELL_MIN_INSTANCE = 0, sigma_gene = FALSE)
saveRDS(RCTD,'../myRCTD_fitted.rds')
RCTD <- readRDS('../myRCTD_fitted.rds')
RCTD@config$max_cores <- 4
testRCTD <- fit.cell.types(RCTD, cell_types)

print(cell.confusion.mat(RCTD,testRCTD)$overall['Accuracy'])


RCTD_results <- run.iter.optim(myRCTD, n_iter = 3)
saveRDS(RCTD_results,'../RCTD_list_TEST.rds')
RCTD_results <- readRDS('../RCTD_list_TEST.rds')

# Run on scRNA-seq and snRNA-seq reference
myRCTD <- readRDS('../cross_platform_RCTD.rds')
RCTD_results <- run.iter.optim(myRCTD, n_iter = 10)
saveRDS(RCTD_results,'../RCTD_list_2.rds')
RCTD_results <- readRDS('../RCTD_list_2.rds')

# Run on DOWNSAMPLED scRNA-seq and snRNA-seq reference
create_downsampled_data <- function(reference, cell_types_keep = NULL, n_samples = 10000) {
	if(is.null(cell_types_keep))
		cell_types_keep = levels(reference@cell_types)
	cell_types_keep = cell_types_keep[unlist(lapply(cell_types_keep, function(x) nchar(x) > 0))]
	index_keep = c(); i = 1
	repeat{
		new_index = which(reference@cell_types == cell_types_keep[i])
		new_samples = min(n_samples, length(new_index))
		index_keep = c(index_keep, sample(new_index,new_samples,replace=FALSE))
		if((i = i + 1) > length(cell_types_keep))
			break
	}
	reference@counts = reference@counts[,index_keep]
	reference@cell_types = reference@cell_types[index_keep]
	reference@cell_types = droplevels(reference@cell_types)
	reference@nUMI = reference@nUMI[index_keep]
	return(reference)
}

train_ref <- readRDS('../reference_RCTD_vec.rds') # reference data
test_ref <- readRDS('../cerebellum_reference_dropviz_spacexr.rds') # data to learn
test_ref <- create_downsampled_data(test_ref, n_samples = 500)
counts <- test_ref@counts
nUMI <- test_ref@nUMI
puck <- SpatialRNA(NULL, counts, nUMI, use_fake_coords = TRUE)
myRCTD <- create.RCTD(puck, train_ref, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
saveRDS(myRCTD,'../cross_platform_RCTD_downsampled.rds')
RCTD_results <- run.iter.optim(myRCTD, n_iter = 10)
saveRDS(RCTD_results,'../RCTD_list_2_downsampled.rds')

myRCTD <- readRDS('../cross_platform_RCTD_downsampled.rds')
RCTD_results <- readRDS('../RCTD_list_2_downsampled.rds')

# ANALYSIS

# Gene type assignments compare to original
for (i in 1:length(RCTD_results)) {
	print(cell.confusion.mat(myRCTD,RCTD_results[[i]]$RCTD)$overall['Accuracy'])
}

# Cell type assignments compare consecutive
for (i in 1:(length(RCTD_results)-1)) {
	print(cell.confusion.mat(RCTD_results[[i]]$RCTD,RCTD_results[[i+1]]$RCTD)$overall['Accuracy'])
}

# Gene expression MSE compare original
for (i in 1:length(RCTD_results)) {
  print(gene.expression.mse(myRCTD,RCTD_results[[i]]$RCTD))
}

# Gene expression MSE compare consecutive
for (i in 1:(length(RCTD_results)-1)) {
  print(gene.expression.mse(RCTD_results[[i]]$RCTD,RCTD_results[[i+1]]$RCTD))
}

# Gene expression correlation heatmap two inputs
correlation <- gene.correlation.mat(RCTD_results[[9]]$RCTD, RCTD_results[[10]]$RCTD)
n_cells <- dim(correlation)[1]/2
correlation <- correlation[1:n_cells, (n_cells+1):(2*n_cells)]
data <- melt(correlation ^ 2)
ggplot(data, aes(Var1, Var2, fill= value)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colors = pals::brewer.blues(20)[2:20], limits=c(0,1), name='r^2') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Gene expression MSE single input
myRCTD <- readRDS('../cross_platform_RCTDE.rds')
gene.expression.mse(myRCTD)

# Gene expression correlation heatmap single input
myRCTD <- readRDS('../myRCTDE.rds')
correlation <- gene.correlation.mat(myRCTD)
n_cells <- dim(correlation)[1]/2
correlation <- correlation[1:n_cells, (n_cells+1):(2*n_cells)]
data <- melt(correlation ^ 2)
ggplot(data, aes(Var1, Var2, fill= value)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colors = pals::brewer.blues(20)[2:20], limits=c(0,1), name='r^2') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



test_ref <- readRDS('../cerebellum_reference_dropviz_spacexr.rds') # data to learn
counts <- test_ref@counts
nUMI <- test_ref@nUMI
puck <- SpatialRNA(NULL, counts, nUMI, use_fake_coords = TRUE)
refRCTD <- create.RCTD(puck, test_ref, max_cores = 4, CELL_MIN_INSTANCE = 0)

# Gene type assignments compare to original reference
for (i in 1:length(RCTD_results)) {
	print(cell.confusion.mat(refRCTD,RCTD_results[[i]]$RCTD)$overall['Accuracy'])
}

# Gene expression MSE compare original
for (i in 1:length(RCTD_results)) {
  print(gene.expression.mse(refRCTD,RCTD_results[[i]]$RCTD))
}

# Find doublet weight statistics
stats_matrix = Matrix(0, nrow = length(RCTD_results), ncol = 3)
colnames(stats_matrix) = c('mean', 'stdev', 'length')
for (i in 1:length(RCTD_results)) {
  stats_matrix[i,] = doublet.weights.statistics(RCTD_results[[i]]$RCTD, spot_types = c('doublet_uncertain'))
}
stats_matrix
