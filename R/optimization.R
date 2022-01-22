library(spacexr)
library(Matrix)
library(doParallel)

#' Refits gene expression profiles for each cell type with CSIDE based on cell type predictions.
#'
#' @param RCTD an RCTD object that has cell types assigned.
#' @param cell_types a vector of cell types whose profiles will be refitted.
#' @param CELL_MIN_INSTANCE minimum number of instances of a cell type for it to be refitted.
#' @param sigma_gene (default TRUE) whether to fit sigma for each gene.
#' @return an RCTD object, with new gene expression profiles stored in de_results.
#' @export
fit.gene.expression <- function(RCTD, cell_types, CELL_MIN_INSTANCE = 25, sigma_gene = TRUE) {
	barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
	X <- rep(1, length(barcodes))
	X <- as.matrix(X)
	rownames(X) <- barcodes
	return(run.CSIDE(RCTD, X, barcodes, cell_types, cell_type_threshold = CELL_MIN_INSTANCE, sigma_gene = sigma_gene,
					  test_genes_sig = FALSE, params_to_test = 1))
}

#' Refits cell type assignments and gene lists with RCTD based on differential expression gene fits.
#'
#' @param RCTD an RCTD object that has differential gene expressions fitted.
#' @param cell_types a vector of cell types to be considered when refitting.
#' @return an RCTD object, with new cell type assignments stored in results.
#' @export
fit.cell.types <- function(RCTD, cell_types) {
	cell_type_info <- list(as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)), cell_types, length(cell_types))
	cell_type_info <- list(info = cell_type_info, renorm = cell_type_info)
	config <- RCTD@config
	puck.original <- RCTD@originalSpatialRNA
	message('fit.cell.types: getting regression differentially expressed genes: ')
	gene_list_reg = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
	if(length(gene_list_reg) == 0)
	  stop("fit.cell.types: Error: 0 regression differentially expressed genes found")
	message('fit.cell.types: getting platform effect normalization differentially expressed genes: ')
	gene_list_bulk = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
	if(length(gene_list_bulk) == 0)
	  stop("fit.cell.types: Error: 0 bulk differentially expressed genes found")
	puck = restrict_counts(puck.original, gene_list_bulk, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
	puck = restrict_puck(puck, colnames(puck@counts))
	RCTD@internal_vars$gene_list_reg <- gene_list_reg; RCTD@internal_vars$gene_list_bulk <- gene_list_bulk
	RCTD@spatialRNA <- puck; RCTD@cell_type_info <- cell_type_info
	return(fitPixels(RCTD))
}

#' Iteratively optimizes cell type assignments and gene expression profiles with RCTD and CSIDE.
#'
#' @param RCTD an RCTD object to run iterative optimization on.
#' @param cell_types_assigned (default FALSE) whether or not RCTD already has cell types assignd in results.
#' @param cell_types (default NULL) cell types to run iterative optimization on. If null, then runs on cell types that appear more than CELL_MIN_INSTANCE times.
#' @param CELL_MIN_INSTANCE (default 25) minimum number of instances of a cell for it to be considered in the optimization.
#' @param max_n_iter (default 20) maximum number of optimization iterations.
#' @param convergence_cutoff (default 0.999) classification accuracy threshold for optimization to terminate early.
#' @param discovery_threshold (default 0.99) value between 0 and 1, threshold of singlet counts between iterations for optimization to terminate early.
#' @param constant_genes (default TRUE) whether or not gene lists change on each iteration.
#' @param used_reference (default FALSE) whether or not a reference was used to generate cell_type_info. If true, the fitBulk is run during the first iteration.
#' @return a list of RCTD objects produced at each iteration of the optimization.
#' @export
run.iter.optim <- function(RCTD, cell_types_assigned = FALSE, cell_types = NULL, CELL_MIN_INSTANCE = 25, max_n_iter = 20, convergence_cutoff = 0.999, discovery_threshold = 0.99, constant_genes = TRUE, used_reference = FALSE) {
	if (!cell_types_assigned) {
		RCTD@config$RCTDmode <- 'doublet'
		if (used_reference)
			RCTD <- fitBulk(RCTD)
		RCTD <- choose_sigma_c(RCTD)
		message('run.iter.optim: assigning initial cell types')
		RCTD <- fitPixels(RCTD, doublet_mode = 'doublet')
	}
	if (is.null(cell_types)) {
		barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
		cell_type_count <- aggregate_cell_types(RCTD, barcodes, doublet_mode = TRUE)
		cell_types <- names(which(cell_type_count >= CELL_MIN_INSTANCE))
	}
	RCTD_list = list()
	RCTD_list[[1]] <- RCTD
	barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
	X <- as.matrix(rep(1, length(barcodes)))
	rownames(X) <- barcodes
	for (i in 1:max_n_iter) {
		message(paste('run.iter.optim: running iteration', i))
		message('fitting gene expression profiles')
		if (!constant_genes) {
			barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
			X <- as.matrix(rep(1, length(barcodes)))
			rownames(X) <- barcodes
		}
		RCTD <- run.CSIDE(RCTD, X, barcodes, cell_types, cell_type_threshold = 0, sigma_gene = FALSE, test_genes_sig = FALSE, params_to_test = 1)
		message('fitting cell types')
		if (constant_genes) {
			new_fits <- as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)); old_fits <- RCTD@cell_type_info$renorm[[1]]
			common_genes <- intersect(rownames(old_fits), rownames(new_fits))
			cell_type_info <- list(rbind(old_fits[! rownames(old_fits) %in% common_genes, cell_types], new_fits[common_genes, cell_types]), cell_types, length(cell_types))
			cell_type_info <- list(info = cell_type_info, renorm = cell_type_info)
			RCTD@cell_type_info <- cell_type_info
		}
		else {
			cell_type_info <- list(as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)), cell_types, length(cell_types))
			cell_type_info <- list(info = cell_type_info, renorm = cell_type_info)
			config <- RCTD@config
			puck.original <- RCTD@originalSpatialRNA
			message('fit.cell.types: getting regression differentially expressed genes: ')
			gene_list_reg = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
			if(length(gene_list_reg) == 0)
			  stop("fit.cell.types: Error: 0 regression differentially expressed genes found")
			message('fit.cell.types: getting platform effect normalization differentially expressed genes: ')
			gene_list_bulk = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
			if(length(gene_list_bulk) == 0)
			  stop("fit.cell.types: Error: 0 bulk differentially expressed genes found")
			puck = restrict_counts(puck.original, gene_list_bulk, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
			puck = restrict_puck(puck, colnames(puck@counts))
			RCTD@internal_vars$gene_list_reg <- gene_list_reg; RCTD@internal_vars$gene_list_bulk <- gene_list_bulk
			RCTD@spatialRNA <- puck; RCTD@cell_type_info <- cell_type_info
		}
		RCTD <- fitPixels(RCTD)
		RCTD_list[[i+1]] <- RCTD
		accuracy <- cell.confusion.mat(RCTD_list[[i]], RCTD_list[[i+1]])$overall['Accuracy']
		n_singlets <- dim(RCTD@results$results_df[RCTD@results$results_df$spot_class == 'singlet',])[1]
		prev_singlets <- dim(RCTD_list[[i]]@results$results_df[RCTD_list[[i]]@results$results_df$spot_class == 'singlet',])[1]
		message(paste('classification convergence accuracy:', accuracy))
		if (accuracy >= convergence_cutoff & discovery_threshold <= n_singlets/prev_singlets & 1/discovery_threshold >= n_singlets/prev_singlets)
			break
	}
	return(RCTD_list)
}
