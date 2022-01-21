library(spacexr)
library(Matrix)
library(doParallel)

#' Creates an RCTD object from a SpatialRNA object without a reference
#'
#' @param spatialRNA a SpatialRNA object to run RCTD on
#' @param gene_cutoff minimum normalized gene expression for genes to be included in the platform effect normalization step.
#' @param fc_cutoff minimum log-fold-change (across cell types) for genes to be included in the platform effect normalization step.
#' @param gene_cutoff_reg minimum normalized gene expression for genes to be included in the RCTD step.
#' @param fc_cutoff_reg minimum log-fold-change (across cell types) for genes to be included in the RCTD step.
#' @param UMI_min minimum UMI per pixel included in the analysis
#' @param UMI_max maximum UMI per pixel included in the analysis
#' @param UMI_min_sigma minimum UMI per pixel for the choose_sigma_c function
#' @param max_cores for parallel processing, the number of cores used. If set to 1, parallel processing is not used. The system will additionally be checked for
#' number of available cores.
#' @param CELL_MIN_INSTANCE minimum number of cells required per cell type. Default 25, can be lowered if desired.
#' @param MAX_MULTI_TYPES (multi-mode only) Default 4, max number of cell types per pixel
#' @param cell_type_info Default NULL, option to pass in cell_type_info directly
#' @param gene_list_reg Default NULL, option to pass in gene_list_reg directly
#' @param gene_list_bulk Default NULL, option to pass in gene_list_bulk directly
#' @param class_df Default NULL, option to pass in class_df directly
#' @return an RCTD object, which is ready to run the run.RCTD function
#' @export
create.RCTD.noref <- function(spatialRNA, max_cores = 4, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 20000000, UMI_min_sigma = 300,
                              CELL_MIN_INSTANCE = 25, MAX_MULTI_TYPES = 4, cell_type_info = NULL, gene_list_reg = NULL, gene_list_bulk = NULL, class_df = NULL) {

	reference <- new("Reference", cell_types = factor(), counts = as(matrix(), 'dgCMatrix'))
	config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma, max_cores = max_cores,
	               N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30, MIN_CHANGE_BULK = 0.0001, MIN_CHANGE_REG = 0.001, UMI_max = UMI_max, MIN_OBS = 3, MAX_MULTI_TYPES = MAX_MULTI_TYPES)
	puck.original = restrict_counts(spatialRNA, rownames(spatialRNA@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
	if(is.null(cell_type_info)) {
		cell_type_info <- list(info = list(data.frame(), c(), 0), renorm = NULL)
		puck = puck.original
	}
	else {
		message('create.RCTD: getting regression differentially expressed genes: ')
		if (is.null(gene_list_reg))
			gene_list_reg = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
		if(length(gene_list_reg) == 0)
		  stop("create.RCTD: Error: 0 regression differentially expressed genes found")
		message('create.RCTD: getting platform effect normalization differentially expressed genes: ')
		if (is.null(gene_list_bulk))
			gene_list_bulk = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
		if(length(gene_list_bulk) == 0)
		  stop("create.RCTD: Error: 0 bulk differentially expressed genes found")
		puck = restrict_counts(puck.original, gene_list_bulk, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
		puck = restrict_puck(puck, colnames(puck@counts))
		if (is.null(class_df))
			class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]]); colnames(class_df)[1] = "class"
	}
	internal_vars <- list(gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk, proportions = NULL, class_df = class_df, cell_types_assigned = F)
	new("RCTD", spatialRNA = puck, originalSpatialRNA = puck.original, reference = reference, config = config, cell_type_info = cell_type_info, internal_vars = internal_vars)
}

#' Sets internal variables of an RCTD object initialized without a reference.
#' 
#' @param RCTD an RCTD object initialized without a reference.
#' @return an RCTD object with internal_vars set.
#' @export
set_internal_vars <- function(RCTD) {
  puck = RCTD@spatialRNA; MIN_UMI = RCTD@config$UMI_min_sigma; sigma = 100
  Q1 <- readRDS(system.file("extdata", "Qmat/Q_mat_1.rds", package = "spacexr"))
  Q2 <- readRDS(system.file("extdata", "Qmat/Q_mat_2.rds", package = "spacexr"))
  Q_mat_all <- c(Q1,Q2)
  sigma_vals <- names(Q_mat_all)
  X_vals <- readRDS(system.file("extdata", "Qmat/X_vals.rds", package = "spacexr"))
  RCTD@internal_vars$sigma <- sigma/100
  RCTD@internal_vars$Q_mat <- Q_mat_all[[as.character(sigma)]]
  RCTD@internal_vars$X_vals <- X_vals
  return(RCTD)
}

#' Creates cell type info from cell type assignments.
#'
#' @param puck a SpatialRNA object to predict cell profiles.
#' @param assignments a dataframe with cell type assignments in column labeled cell_types and rows labeled with barcodes.
#' @return a list of cell type info that can be inputted into create.RCTD.
#' @export
cell_type_info_from_assignments <- function(puck, assignments) {
	assigned_cell_types <- assignments$cell_types
	names(assigned_cell_types) <- rownames(assignments)
	info <- get_cell_type_info(puck@counts, assigned_cell_types, puck@nUMI)
	return(list(info = info, renorm = info))
}

#' Creates uniformly random generated cell type info.
#'
#' @param puck a SpatialRNA object to predict cell profiles.
#' @param n_cell_types the number of cell types to be assumed in cell type info.
#' @return a list of cell type info that can be inputted into create.RCTD.
#' @export
gen.random.info <- function(puck, n_cell_types) {
	genes <- rownames(puck@counts)
	info <- data.frame(replicate(n_cell_types,runif(length(genes))))
	info <- apply(info, 2, function(x) x/sum(x))
	colnames(info) <- 1:n_cell_types; rownames(info) <- genes
	return(list(info = list(info, 1:n_cell_types, n_cell_types), renorm = list(info, 1:n_cell_types, n_cell_types)))
}

#' Records cell types in RCTD object results from cell type assignments.
#'
#' @param RCTD an RCTD object where cell types will be assigned.
#' @param assignments a dataframe with cell type assignments in column labeled cell_types and rows labeled with barcodes.
#' @param weight (default 0.9) the weight given to the assigned cell type in the weights_doublet matrix.
#' @return an RCTD object with cell type predictions in results.
#' @export
assign.cell.types <- function(RCTD, assignments, weight = 0.9) {
	RCTD@internal_vars$cell_types_assigned <- TRUE
	RCTD@config$RCTDmode <- "doublet"
	cell_type_names = levels(assignments[,1])
	RCTD@cell_type_info <- list(info = list(data.frame(), cell_type_names, length(cell_type_names)), renorm = NULL)
	barcodes <- colnames(RCTD@spatialRNA@counts)
	N <- length(barcodes)
	weights_doublet = Matrix(0, nrow = N, ncol = 2)
	rownames(weights_doublet) = barcodes
	colnames(weights_doublet) = c('first_type', 'second_type')
	empty_cell_types = factor(character(N),levels = cell_type_names)
	spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_uncertain")
	results_df <- data.frame(spot_class = factor(character(N),levels=spot_levels), first_type = empty_cell_types, second_type = empty_cell_types)
	for(i in 1:N) {
		if(i %% 1000 == 0)
			print(paste("assign.cell.types: finished",i))
	  barcode <- barcodes[i]
		weights_doublet[i,] = c(weight, 1-weight)
		results_df[i, "spot_class"] = "singlet"
		results_df[i, "first_type"] = assignments[barcode,1]
		results_df[i, "second_type"] = assignments[barcode,1]
	}
	rownames(results_df) = barcodes
	RCTD@results <- list(results_df = results_df, weights_doublet = weights_doublet)
  	return(RCTD)
}

