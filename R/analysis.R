library(spacexr)
library(Matrix)
library(caret)
library(reshape2)

#' Plots spatial location of predicted cell types and stores results in UROP/results
#'
#' @param RCTD an RCTD object that has cell types stored in results.
#' @return a plot of the spatial locations of cell types
#' @export
plot.cell.types <- function(RCTD) {
	plot_all_cell_types(RCTD@results$results_df, RCTD@originalSpatialRNA@coords, RCTD@cell_type_info$renorm[[2]], '..')
}

#' Compares singlet cell type assignments by generating a confusion matrix using the caret package.
#' 
#' @param RCTD_ref an RCTD object that has ground truth cell types stored in results.
#' @param RCTD_pred an RCTD object that has predicted cell types stored in results.
#' @param cell_types (default NULL) a vector of cell types included in the confusion matrix. If null, cell types will be all predicted cell types.
#' @param verbose (default FALSE) whether to display number of reference, predicted, and common singlets.
#' @return a list object that contains the confusion matrix and related statistics.
#' @export
cell.confusion.mat <- function(RCTD_ref, RCTD_pred, cell_types = NULL, verbose = FALSE) {
	if (is.null(cell_types))
		cell_types <- RCTD_pred@internal_vars_de$cell_types_present
	ref_singlet <- RCTD_ref@results$results_df
	ref_singlet <- ref_singlet[ref_singlet$spot_class == 'singlet' & is.element(ref_singlet$first_type, cell_types),]
	pred_singlet <- RCTD_pred@results$results_df
	pred_singlet <- pred_singlet[pred_singlet$spot_class == 'singlet' & is.element(pred_singlet$first_type, cell_types),]
	common_barcode <- intersect(row.names(ref_singlet), row.names(pred_singlet))
	if (verbose) {
		message(paste(dim(ref_singlet)[1], 'reference singlets'))
		message(paste(dim(pred_singlet)[1], 'predicted singlets'))
		message(paste(length(common_barcode), 'common singlets'))
	}
	ref_singlet <- ref_singlet[common_barcode,]
	pred_singlet <- pred_singlet[common_barcode,]
	ref_types = unlist(list(ref_singlet[,'first_type']))
	pred_types = unlist(list(pred_singlet[,'first_type']))
	return(confusionMatrix(pred_types, factor(ref_types, levels = levels(pred_types))))
}

#' Compares singlet cell type assignments by generating a two way table of cell types.
#' 
#' @param RCTD_truth an RCTD object that has ground truth cell types stored in results.
#' @param RCTD_pred an RCTD object that has predicted cell types stored in results.
#' @param verbose (default FALSE) whether to display number of reference, predicted, and common singlets.
#' @return a table containing cell type assignments.
#' @export
cell.table <- function(RCTD_truth, RCTD_pred, verbose = FALSE) {
  truth_singlet <- RCTD_truth@results$results_df
  truth_singlet <- truth_singlet[truth_singlet$spot_class == 'singlet',]
  pred_singlet <- RCTD_pred@results$results_df
  pred_singlet <- pred_singlet[pred_singlet$spot_class == 'singlet',]
  common_barcode <- intersect(row.names(truth_singlet), row.names(pred_singlet))
  if (verbose) {
		message(paste(dim(truth_singlet)[1], 'ground truth singlets'))
	  message(paste(dim(pred_singlet)[1], 'predicted singlets'))
	  message(paste(length(common_barcode), 'common singlets'))
  }
  truth_singlet <- truth_singlet[common_barcode,]
  pred_singlet <- pred_singlet[common_barcode,]
  truth_types = unlist(list(truth_singlet[,'first_type']))
  pred_types = unlist(list(pred_singlet[,'first_type']))
  return(table(truth_types, pred_types))
}

#' Generates a heatmap plot of cell type classification accuracy between two RCTD objects.
#'
#' @param RCTD_truth an RCTD object that has ground truth cell types stored in results.
#' @param RCTD_pred an RCTD object that has predicted cell types stored in results.
#' @param cell_table (default NULL) option to pass in cell assignment table directly
#' @param verbose (default FALSE) whether to display number of reference, predicted, and common singlets.
#' @param cell_min_instance (default 10) minimum number of ground truth cells required to be included in the heatmap.
#' @return a heatmap plot of the cell type classification accuracy.
#' @export
cell.heatmap <- function(RCTD_truth, RCTD_pred, mytable = NULL, verbose = FALSE, cell_min_instance = 10) {
	if (is.null(mytable))
		mytable <- cell.table(RCTD_truth, RCTD_pred, verbose = verbose)
	mytable <- mytable[rowSums(mytable) >= cell_min_instance, ]
	mytable <- apply(mytable, 1, function(x) x/sum(x))
	data <- melt(mytable)
	plot <- ggplot(data, aes(factor(truth_types), factor(pred_types), fill= value)) +
							   geom_tile() +
							   theme_classic() +
							   scale_fill_gradientn(colors = pals::brewer.blues(20)[2:20], limits=c(0,1), name='Classification Proportion') +
							   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
							   xlab('Reference Cell Type')+ ylab('Predicted Cell Type')
	return(plot)
}

#' Calculates conditional entropy between cell type assignments in two different RCTD objects.
#'
#' The conditional entropy of Y given X is denoted as H(Y|X). Computes H(RCTD_pred|RCTD_truth).
#' 
#' @param RCTD_truth an RCTD object that has ground truth cell types stored in results.
#' @param RCTD_pred an RCTD object that has predicted cell types stored in results.
#' @param cell_table (default NULL) option to pass in cell assignment table directly
#' @param verbose (default FALSE) whether to display number of reference, predicted, and common singlets.
#' @return a value corresponding to conditional entropy.
#' @export
cond.entropy <- function(RCTD_truth, RCTD_pred, mytable = NULL, verbose = FALSE) {
	if (is.null(mytable))
		mytable <- cell.table(RCTD_truth, RCTD_pred, verbose = verbose)
	mytable <- mytable/sum(mytable) * log(mytable/rowSums(mytable))
	mytable[is.nan(mytable)] <- 0
	return(-sum(mytable))
}

#' Calculates classification accuracy between cell type assignments in two different RCTD objects.
#' 
#' @param RCTD_truth an RCTD object that has ground truth cell types stored in results.
#' @param RCTD_pred an RCTD object that has predicted cell types stored in results.
#' @param cell_table (default NULL) option to pass in cell assignment table directly
#' @param verbose (default FALSE) whether to display number of reference, predicted, and common singlets.
#' @return a value corresponding to classification accuracy.
#' @export
class.accuracy <- function(RCTD_truth, RCTD_pred, mytable = NULL, verbose = FALSE) {
	if (is.null(mytable))
		mytable <- cell.table(RCTD_truth, RCTD_pred, verbose = verbose)
	return(sum(apply(mytable, 1, max))/sum(mytable))
}

#' Generates plot of classification of accuracy by ground truth cell type.
#' 
#' @param RCTD_truth an RCTD object that has ground truth cell types stored in results.
#' @param RCTD_pred an RCTD object that has predicted cell types stored in results.
#' @param cell_table (default NULL) option to pass in cell assignment table directly
#' @param verbose (default FALSE) whether to display number of reference, predicted, and common singlets.
#' @param cell_min_instance (default 10) minimum number of ground truth cells required to be included in the heatmap.
#' @return a plot of classification accuracy vs. cell type.
#' @export
class.accuracy.plot <- function(RCTD_truth, RCTD_pred, mytable = NULL, verbose = FALSE, cell_min_instance = 10) {
	if (is.null(mytable))
		mytable <- cell.table(RCTD_truth, RCTD_pred, verbose = verbose)
	mytable <- mytable[rowSums(mytable) >= cell_min_instance, ]
	mydf <- as.data.frame(apply(mytable, 1, max) / rowSums(mytable))
	colnames(mydf) <- 'proportion'
	mydf$cell_type <- factor(rownames(mydf))
	p <- ggplot(mydf, aes(x=cell_type, y=proportion)) + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('Ground Truth Cell Type')+ ylab('Classification Proportion') +
    ylim(0,1)
  return(p)
}

assignment.accuracy <- function(mytable, N = NULL) {
	cluster_assignments <- apply(mytable, 2, function(x) which(x == max(x)))
	truth_totals <- c(dim(mytable)[1])
	for (i in 1:dim(mytable)[1]) {
		truth_totals[i] <- sum(mytable[i,] * (cluster_assignments == i))
	}
	cell_accuracy <- truth_totals / apply(mytable, 1, sum)
	cell_accuracy <- sort(cell_accuracy, decreasing=TRUE)
	cell_accuracy <- cell_accuracy[cell_accuracy != 0]
	if (is.null(N))
		N = length(cell_accuracy)
	else 
		N = min(length(cell_accuracy), N)
	return(mean(cell_accuracy[1:N]))
}

assignment.accuracy.plot <- function(mytable, cell_min_instance = 1) {
	mytable <- mytable[rowSums(mytable) >= cell_min_instance, ]
	cluster_assignments <- apply(mytable, 2, function(x) which(x == max(x)))
	truth_totals <- c(dim(mytable)[1])
	for (i in 1:dim(mytable)[1]) {
		truth_totals[i] <- sum(mytable[i,] * (cluster_assignments == i))
	}
	mydf <- as.data.frame(truth_totals / apply(mytable, 1, sum))
	colnames(mydf) <- 'proportion'
	mydf$cell_type <- factor(rownames(mydf))
	mydf <- mydf[mydf$proportion != 0,]
	p <- ggplot(mydf, aes(x=cell_type, y=proportion)) + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('Ground Truth Cell Type')+ ylab('Classification Proportion') +
    ylim(0,1)
  return(p)
}

#' Finds gene expression profiles of two RCTD objects.
#'
#' @param RCTD_ref an RCTD object with reference gene expression stored in cell_type_info.
#' @param RCTD_pred an RCTD object with predicted gene expression stored in cell_type_info.
#' @param gene_list (default NULL) a vector of genes to be used when computing gene expression profiles. If null, selects all genes shared between the two RCTD objects.
#' @return a list of the reference gene expression profile and predicted gene expression profile.
#' @export
gene.expression <- function(RCTD_ref, RCTD_pred, gene_list = NULL) {
	pred_expression <- RCTD_pred@cell_type_info$renorm[[1]]
	ref_expression <- RCTD_ref@cell_type_info$renorm[[1]]
	if (is.null(gene_list))
		gene_list <- intersect(rownames(na.omit(ref_expression)), rownames(na.omit(pred_expression)))
	ref_expression <- data.matrix(ref_expression[gene_list,])
	pred_expression <- data.matrix(pred_expression[gene_list,])
	return(list(reference = ref_expression, prediction = pred_expression))
}

#' Generates gene expression correlation matrix of log gene expression between two RCTD objects.
#'
#' @param RCTD_ref an RCTD object with reference gene expression stored in cell_type_info.
#' @param RCTD_pred an RCTD object with predicted gene expression stored in cell_type_info.
#' @param gene_list (default NULL) a vector of genes to be used when computing gene expression profiles. If null, selects all genes shared between the two RCTD objects.
#' @param relabel (default FALSE) whether or not to rename cell types to distinguish between reference and predicted.
#' @return a matrix of the correlation between gene expression profiles between cell types.
#' @export
gene.correlation.mat <- function(RCTD_ref, RCTD_pred, gene_list = NULL, relabel = FALSE) {
	gene_expression <- gene.expression(RCTD_ref, RCTD_pred, gene_list = gene_list)
	ref_expression <- log(gene_expression$reference); pred_expression <- log(gene_expression$prediction)
	ref_cells <- colnames(ref_expression); pred_cells <- colnames(pred_expression)
	n_ref_cells <- length(ref_cells); n_pred_cells <- length(pred_cells)
	if (relabel) {
		colnames(ref_expression) <- lapply(ref_cells, function(x) paste0(x, '_ref'))
		colnames(pred_expression) <- lapply(pred_cells, function(x) paste0(x, '_pred'))
	}
	expression <- as.matrix(cbind(ref_expression, pred_expression))
	expression <- expression[rowSums(is.infinite(expression)) == 0, ]
	return(cor(expression)[1:n_ref_cells,(n_ref_cells + 1):(n_ref_cells + n_pred_cells)])
}

#' Generates a heatmap plot of gene expression correlation by cell type between two RCTD objects.
#'
#' @param RCTD_ref an RCTD object with reference gene expression stored in cell_type_info.
#' @param RCTD_pred an RCTD object with predicted gene expression stored in cell_type_info.
#' @param correlation (default NULL) option to pass in correlation matrix directly
#' @param gene_list (default NULL) a vector of genes to be used when computing gene expression profiles. If null, selects all genes shared between the two RCTD objects.
#' @return a heatmap plot of the correlation between gene expression profiles between cell types.
#' @export
gene.heatmap <- function(RCTD_ref, RCTD_pred, correlation = NULL, gene_list = NULL) {
	if (is.null(correlation))
		correlation <- gene.correlation.mat(RCTD_ref, RCTD_pred, gene_list = gene_list)
	data <- melt(correlation ^ 2)
	plot <- ggplot(data, aes(factor(Var1), factor(Var2), fill= value)) +
							   geom_tile() +
							   theme_classic() +
							   scale_fill_gradientn(colors = pals::brewer.blues(20)[2:20], limits=c(0,1), name='r^2') +
							   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
							   xlab('Reference Cell Type')+ ylab('Predicted Cell Type')
	return(plot)
}

#' Computes mean squared error between log gene expression profiles using cell types in common.
#'
#' @param RCTD_ref an RCTD object with reference gene expression stored in cell_type_info.
#' @param RCTD_pred an RCTD object with predicted gene expression stored in cell_type_info.
#' @param gene_list (default NULL) a vector of genes to be used when computing gene expression profiles. If null, selects all genes shared between the two RCTD objects.
#' @return a value corresponding to the mean squared error between log gene expression profiles.
#' @export
gene.mse <- function(RCTD_ref, RCTD_pred, gene_list = NULL) {
	pred_expression <- data.matrix(log(RCTD_pred@cell_type_info$renorm[[1]]))
	ref_expression <- data.matrix(log(RCTD_ref@cell_type_info$renorm[[1]]))
	common_cells <- intersect(colnames(ref_expression), colnames(pred_expression))
	if (is.null(gene_list))
		gene_list <- intersect(rownames(na.omit(ref_expression[rowSums(is.infinite(ref_expression)) == 0,])), rownames(na.omit(pred_expression[rowSums(is.infinite(pred_expression)) == 0,])))
	ref_expression <- ref_expression[gene_list,common_cells]
	pred_expression <- pred_expression[gene_list,common_cells]
	return(sum((ref_expression - pred_expression)^2) / length(ref_expression))
}

#' Computes mean squared error between log gene expression profiles for each pair of cell types in two RCTD objects.
#'
#' @param RCTD_ref an RCTD object with reference gene expression stored in cell_type_info.
#' @param RCTD_pred an RCTD object with predicted gene expression stored in cell_type_info.
#' @param gene_list (default NULL) a vector of genes to be used when computing gene expression profiles. If null, selects all genes shared between the two RCTD objects.
#' @return a matrix of mean squared error between log gene expression profiles between cell types.
#' @export
gene.mse.mat <- function(RCTD_ref, RCTD_pred, gene_list = NULL) {
	pred_expression <- data.matrix(log(RCTD_pred@cell_type_info$renorm[[1]]))
	ref_expression <- data.matrix(log(RCTD_ref@cell_type_info$renorm[[1]]))
	if (is.null(gene_list))
		gene_list <- intersect(rownames(na.omit(ref_expression[rowSums(is.infinite(ref_expression)) == 0,])), rownames(na.omit(pred_expression[rowSums(is.infinite(pred_expression)) == 0,])))
	ref_expression <- ref_expression[gene_list,]
	pred_expression <- pred_expression[gene_list,]
	ref_cells <- colnames(ref_expression); pred_cells <- colnames(pred_expression)
	mse_mat <- matrix(,nrow = length(ref_cells), ncol = length(pred_cells))
	rownames(mse_mat) <- ref_cells; colnames(mse_mat) <- pred_cells
	for (i in 1:length(ref_cells)) {
		for (j in 1:length(pred_cells)) {
			mse_mat[i,j] <- sum((ref_expression[,i] - pred_expression[,j])^2) / length(gene_list)
		}
	}
	return(mse_mat)
}

gene.log.mse <- function(gene_expression, truth_cell, pred_cell) {
	ref_expression <- log(gene_expression$reference[,truth_cell])
	pred_expression <- log(gene_expression$prediction[,pred_cell])
	diff_expression <- ref_expression - pred_expression
	diff_expression <- na.omit(diff_expression[is.finite(diff_expression)])
	return(sum(diff_expression^2)/ length(diff_expression))
}

gene.aggregate.mse <- function(gene_expression, cell_table) {
	cluster_assignments <- apply(cell_table, 2, function(x) which(x == max(x)))
	tot_mse = 0
	for (i in 1:length(cluster_assignments)) {
		pred_cell = i
		truth_cell = cluster_assignments[i]
		mse <- gene.log.mse(gene_expression, truth_cell, pred_cell)
		tot_mse = tot_mse + mse
	}
	return(tot_mse/length(cluster_assignments))
}

gene.expression.rank <- function(gene_expression, truth_cell, pred_cell) {
	ref_expression <- sort(gene_expression$reference[,truth_cell], decreasing = TRUE)
	pred_expression <- sort(gene_expression$prediction[,pred_cell], decreasing = TRUE)
	for (i in 1:length(ref_expression)) {
	    ref_expression[i] = i
	}
	for (i in 1:length(pred_expression)) {
	    pred_expression[i] = i
	}
	ref_expression <- ref_expression[order(factor(names(ref_expression)))]
	pred_expression <- pred_expression[order(factor(names(pred_expression)))]
	return(list(reference = ref_expression, prediction = pred_expression))
}

aggregate.gene.rank <- function(gene_expression, cell_table) {
	cluster_assignments <- apply(cell_table, 2, function(x) which(x == max(x)))
	pred_ranks = c()
	ref_ranks = c()
	for (i in 1:length(cluster_assignments)) {
		pred_cell = i
		truth_cell = cluster_assignments[i]
		rank <- gene.expression.rank(gene_expression, truth_cell, pred_cell)
		pred_ranks <- append(pred_ranks, rank$prediction)
		ref_ranks <- append(ref_ranks, rank$reference)
	}
	return(list(reference = ref_ranks, prediction = pred_ranks))
}

gene.bin.plot <- function(rank, rank_lim=NA, bins = 20) {
	rank <- as.data.frame(rank)
	p <- ggplot(rank, aes(x=reference, y=prediction)) +
	  geom_bin2d(bins=bins) +
	  scale_fill_continuous(type = "viridis") +
	  theme_bw() +
	  xlab('Reference Gene Expression Rank') + ylab('Predicted Gene Expression Rank') +
	  xlim(c(NA, rank_lim)) + ylim(c(NA, rank_lim))
	return(p)
}

gene.density.plot <- function(rank, rank_lim=NA) {
	rank <- as.data.frame(rank)
	p <- ggplot(rank, aes(x=reference, y=prediction) ) +
	  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
	  scale_fill_distiller(palette= "Spectral", direction=1) +
	  scale_x_continuous(expand = c(0, 0)) +
	  scale_y_continuous(expand = c(0, 0)) +
	  xlab('Reference Gene Expression Rank') + ylab('Predicted Gene Expression Rank') +
	  xlim(c(NA, rank_lim)) + ylim(c(NA, rank_lim))
	return(p)
}

#' Computes statistics of doublet weights.
#'
#' @param RCTD an RCTD object that has cell types predicted.
#' @param spot_types a vector of spot types to find statistics for.
#' @return a vector of the mean, standard deviation, and number of spots.
#' @export
doublet.weights.statistics <- function(RCTD, spot_types) {
	results_df <- RCTD@results$results_df
	doublet_weights <- RCTD@results$weights_doublet
	prob <- c()
	for (i in 1:dim(results_df)[1]) {
	    if (results_df[i, 'spot_class'] %in% spot_types)
	        prob <- append(prob, doublet_weights[rownames(results_df)[i],1])
	}
	return(c(mean(prob), sd(prob), length(prob)))
}
