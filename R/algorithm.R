library(caret)
library(doParallel)
library(Matrix)
library(Seurat)
library(spacexr)
library(cluster)


#' Unsupervised assignment of cell types.
#'
#' Performs unsupervised cell type assignment and gene expression prediction
#' for spatial transcriptomics data using an iterative optimization algorithm
#' with a clustering-based initialization.
#'
#' @param puck a SpatialRNA object to run prediction on
#' @param max_cores (default 4) maximum number of cores for parallel processing
#' @param max_iter (default 50) maximum number of optimization iterations
#' @param resolution (default 1) resolution for initial clustering
#' @param doublet_mode (default 'doublet') mode to run the prediction
#'   ('doublet' or 'full')
#' @param use_silhouette (default FALSE) whether to use silhouette scores to
#'   trim the initial clustering
#' @param silhouette_cutoff fraction of clustered spots to trim
#' @param SCT (default TRUE) whether to use SCTransform for initial clustering
#' @param gene_list (default NULL) option to input a custom gene list to run 
#'   the prediction. By default, calculates differentially-expressed genes
#' @param info_type (default 'mean') method of initial gene expression profile
#'   generation ('mean' or 'de')
#' @param fit_genes (default 'de') method of gene expression profile generation
#'   ('mean', 'de', or 'singlet de')
run.unsupervised <- function(puck,
  													 max_cores = 4,
  													 max_iter = 50,
  													 resolution = 1,
  													 doublet_mode = 'doublet',
  													 use_silhouette = F,
  													 silhouette_cutoff = 0,
  													 SCT = T,
  													 gene_list = NULL,
  													 info_type = 'mean',
  													 fit_genes = 'de') {
  assignments <- gen.clusters(
  	puck,
  	resolution = resolution,
  	use_silhouette = use_silhouette,
  	silhouette_cutoff = silhouette_cutoff,
  	SCT = SCT
  )
  if (info_type == 'mean') {
    cell_type_info <- cell_type_info_from_assignments(puck, assignments)
  } else if (info_type == 'de') {
    cell_type_info <- de_info_from_assignments(puck, assignments)
  } else {
    stop(paste0(
    	"run.unsupervised: info_type=",
    	info_type,
    	" is not a valid choice. Please set info_type=mean or de."))
  }
  myRCTD <- create.RCTD.noref(
  	puck,
  	max_cores = max_cores,
  	cell_type_info = cell_type_info,
  	gene_list_reg = gene_list
  )
  iter.optim(
  	myRCTD,
  	doublet_mode = doublet_mode,
  	fit_genes = fit_genes,
  	max_iter = max_iter
  )
}


#' Semisupervised assignment of cell types.
#'
#' Performs semisupervised cell type assignment and gene expression prediction
#' for spatial transcriptomics data using an iterative optimization algorithm
#' provided a given initialization.
#'
#' @param RCTD a labeled RCTD object to run prediction on
#' @param max_cores (default 4) maximum number of cores for parallel processing
#' @param max_iter (default 50) maximum number of optimization iterations
#' @param doublet_mode (default 'doublet') mode to run the prediction
#'   ('doublet' or 'full')
#' @param gene_list (default NULL) option to input a custom gene list to run 
#'   the prediction. By default, calculates highly-expressed genes
#' @param fit_genes (default 'de') method of gene expression profile generation
#'   ('mean', 'de', or 'singlet de')
run.semisupervised <- function(RCTD,
														   max_cores = 4, 
														   max_iter = 50,
														   doublet_mode = 'doublet',
														   gene_list = NULL,
														   fit_genes = 'de') {
	if (is.null(gene_list)) {
		gene_threshold = 5e-5
	} else {
		RCTD@originalSpatialRNA <- restrict_counts(RCTD@originalSpatialRNA, 
			gene_list,
			UMI_thresh = 100, 
			UMI_max = 20000000
		)
		gene_threshold = 0
	}
	cell_type_info <- fit.gene.expression(RCTD, gene_threshold = gene_threshold)
	myRCTD <- create.RCTD.noref(
		RCTD@originalSpatialRNA,
		max_cores = max_cores,
		cell_type_info = cell_type_info, 
		gene_list_reg = gene_list,
		class_df = RCTD@internal_vars$class_df
	)
	iter.optim(
		myRCTD,
		doublet_mode = doublet_mode,
		fit_genes = fit_genes,
		max_iter = max_iter,
	)
}


iter.optim <- function(RCTD, 
											 doublet_mode = 'doublet',
											 fit_genes = 'de',
											 cell_types = NULL, 
											 max_iter = 50,
											 convergence_thresh = 0.999) {
	RCTD@config$RCTDmode <- doublet_mode
	RCTD <- choose_sigma_c(RCTD)
	message('iter.optim: assigning initial cell types')
	RCTD <- fitPixels(RCTD, doublet_mode = doublet_mode)
	barcodes <- intersect(
		names(RCTD@spatialRNA@nUMI), 
		colnames(RCTD@spatialRNA@counts)
	)
	if (is.null(cell_types)) cell_types <- RCTD@cell_type_info$info[[2]]
	RCTD_prev <- RCTD
	X <- as.matrix(rep(1, length(barcodes)))
	rownames(X) <- barcodes
	for (i in 1:max_iter) {
		message(paste('iter.optim: running iteration', i))
		message('fitting gene expression profiles')
		if (doublet_mode == 'doublet') {
			if (fit_genes == 'de') {
	  		RCTD <- run.CSIDE(
	  			RCTD,
	  			X,
	  			barcodes,
	  			cell_types,
	  			cell_type_threshold = 0,
	  			gene_threshold = 0,
	  			sigma_gene = F,
	  			test_genes_sig = F,
	  			params_to_test = 1
	  		)
	  		info <- list(
	  			as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)),
	  			cell_types,
	  			length(cell_types)
	  		)
			} else if (fit_genes == 'mean') {
				results <- RCTD@results$results_df
			  pred_singlets <- results[results$spot_class == 'singlet',]
			  cell_type_list <- pred_singlets$first_type
			  names(cell_type_list) <- rownames(pred_singlets)
			  info <- get_cell_type_info(
			  	RCTD@originalSpatialRNA@counts[, rownames(pred_singlets)],
			  	cell_type_list,
			  	RCTD@originalSpatialRNA@nUMI[rownames(pred_singlets)]
			  )
			} else if (fit_genes == 'singlet de') {
				results <- RCTD@results$results_df
				results[results$spot_class != 'singlet', ]$spot_class <- 'reject'
				RCTD <- run.CSIDE(
					RCTD,
					X,
					barcodes,
					cell_types,
					cell_type_threshold = 0,
					gene_threshold = 0,
					sigma_gene = F,
					test_genes_sig = F,
					params_to_test = 1
				)
	  		info <- list(
	  			as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)),
	  			cell_types,
	  			length(cell_types)
	  		)
			} else {
			  stop(paste0(
			  	"iter.optim: fit_genes=", 
			  	fit_genes,
			  	" is not a valid choice. ",
			  	"Please set fit_genes=mean, de, or singlet de.")
			  )
			}
		} else if (doublet_mode == 'full') {
			RCTD <- run.CSIDE(
				RCTD,
				X,
				barcodes,
				cell_types,
				doublet_mode = F,
				cell_type_threshold = 0,
				gene_threshold = 0,
				sigma_gene = F,
				test_genes_sig = F,
				params_to_test = 1
			)
  		info <- list(
  			as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)),
  			cell_types,
  			length(cell_types)
  		)
		}
		RCTD@cell_type_info <- list(info = info, renorm = info)
		message('fitting cell types')
		RCTD <- fitPixels(RCTD, doublet_mode = doublet_mode)
		if (doublet_mode == 'doublet') {
			accuracy <- assignment_accuracy(RCTD_prev, RCTD)
			message(paste('pixel assignment accuracy:', accuracy))
			if (accuracy > convergence_thresh) break
		}
		if (doublet_mode == 'full') {
			mse <- weights_mse(RCTD_prev, RCTD)
			message(paste('pixel weights mse:', mse))
			if (mse < 1e-6) break
		}
		RCTD_prev <- RCTD
	}
	if (doublet_mode == 'doublet') {
		results <- RCTD@results$results_df
		reassign <- rownames(
			results[results$spot_class == 'doublet_certain' & 
			results$singlet_score - results$min_score < 25, ]
		)
		RCTD@results$results_df[reassign, ]$spot_class <- 'singlet'
	}
	RCTD
}


assignment_accuracy <- function(RCTD1, RCTD2) {
  results1 <- RCTD1@results$results_df[, 1:3]
  results2 <- RCTD2@results$results_df[, 1:3]
  placeholder <- levels(results1$first_type)[1]
  results1[results1$spot_class == 'singlet', ]$second_type <- placeholder
  results1[results1$spot_class == 'reject', ]$first_type <- placeholder
  results1[results1$spot_class == 'reject', ]$second_type <- placeholder
  results2[results2$spot_class == 'singlet', ]$second_type <- placeholder
  results2[results2$spot_class == 'reject', ]$first_type <- placeholder
  results2[results2$spot_class == 'reject', ]$second_type <- placeholder
  common <- as.data.frame(results1 == results2)
  correct <- dim(
  	common[common$spot_class &
  	common$first_type &
  	common$second_type, ])[1]
  correct / dim(results1)[1]
}


weights_mse <- function(RCTD1, RCTD2) {
	weights1 <- RCTD1@results$weights
	weights2 <- RCTD2@results$weights
	sum((log(weights1) - log(weights2)) ** 2) / length(weights1)
}


create.RCTD.noref <- function(spatialRNA,
															max_cores = 4,
															gene_cutoff_reg = 0.0002,
															fc_cutoff_reg = 0.75,
															UMI_min = 100,
															UMI_max = 20000000, 
							  							UMI_min_sigma = 300,
							  							MAX_MULTI_TYPES = 4,
							  							cell_type_info = NULL,
							  							gene_list_reg = NULL,
							  							class_df = NULL) {
	config <- list(
		gene_cutoff_reg = gene_cutoff_reg,
		fc_cutoff_reg = fc_cutoff_reg,
		UMI_min = UMI_min,
		UMI_min_sigma = UMI_min_sigma,
		max_cores = max_cores,
		N_epoch = 8,
		N_X = 50000,
		K_val = 100,
		N_fit = 1000,
		N_epoch_bulk = 30,
		MIN_CHANGE_BULK = 0.0001,
		MIN_CHANGE_REG = 0.001,
		UMI_max = UMI_max,
		MIN_OBS = 3,
		MAX_MULTI_TYPES = MAX_MULTI_TYPES
	)
	reference <- new(
		"Reference",
		cell_types = factor(),
		counts = as(matrix(), 'dgCMatrix')
	)
	puck.original = restrict_counts(spatialRNA,
		rownames(spatialRNA@counts),
		UMI_thresh = config$UMI_min,
		UMI_max = config$UMI_max
	)
	if (is.null(cell_type_info)) {
		cell_type_info <- list(
			info = list(data.frame(), c(), 0),
			renorm = list(data.frame(), c(), 0)
		)
		puck = puck.original
	} else {
		if (is.null(gene_list_reg)) {
			message(paste(
				"create.RCTD.noref:",
				"getting regression differentially expressed genes: "
			))
			gene_list_reg = get_de_genes(
				cell_type_info$info,
				puck.original,
				fc_thresh = config$fc_cutoff_reg,
				expr_thresh = config$gene_cutoff_reg,
				MIN_OBS = config$MIN_OBS
			)
		}
		if(length(gene_list_reg) == 0) {
		  stop(paste(
		  	"create.RCTD.noref: Error:",
		  	"0 regression differentially expressed genes found"
		  ))
		}
		puck = restrict_counts(
			spatialRNA,
			gene_list_reg,
			UMI_thresh = config$UMI_min,
			UMI_max = config$UMI_max
		)
		puck = restrict_puck(puck, colnames(puck@counts))
		puck.original = puck
		if (is.null(class_df))
			class_df <- data.frame(
				cell_type_info$info[[2]],
				row.names = cell_type_info$info[[2]]
			)
			colnames(class_df)[1] = "class"
		cell_type_info$info[[1]] <- cell_type_info$info[[1]][gene_list_reg, ]
		cell_type_info$renorm[[1]] <- cell_type_info$renorm[[1]][gene_list_reg, ]
	}
	internal_vars <- list(
		gene_list_reg = gene_list_reg,
		proportions = NULL,
		class_df = class_df,
		cell_types_assigned = F
	)
	new(
		"RCTD",
		spatialRNA = puck,
		originalSpatialRNA = puck.original,
		reference = reference,
		config = config,
		cell_type_info = cell_type_info,
		internal_vars = internal_vars
	)
}


cell_type_info_from_assignments <- function(puck, assignments) {
	counts <- puck@counts[, rownames(assignments)]
	nUMI <- puck@nUMI[rownames(assignments)]
	assigned_cell_types <- assignments$cell_types
	names(assigned_cell_types) <- rownames(assignments)
	info <- get_cell_type_info(counts, assigned_cell_types, nUMI)
	list(info = info, renorm = info)
}


de_info_from_assignments <- function(puck,
																		 assignments,
																		 gene_threshold = 5e-5) {
	myRCTD <- create.RCTD.noref(puck)
	myRCTD <- set_internal_vars(myRCTD)
	myRCTD <- assign.cell.types(myRCTD, assignments)
	fit.gene.expression(myRCTD, gene_threshold = gene_threshold)
}


gen.clusters <- function(puck,
												 resolution = 1,
												 SCT = T,
												 use_silhouette = F,
												 silhouette_cutoff = 0) {
	slide.seq <- CreateSeuratObject(counts = puck@counts, assay = "Spatial")
	if (SCT) {
		slide.seq <- SCTransform(slide.seq, assay = "Spatial", verbose = F)
		slide.seq <- RunPCA(slide.seq, assay = "SCT")
	} else {
		slide.seq <- NormalizeData(slide.seq)
		slide.seq <- FindVariableFeatures(slide.seq)
		slide.seq <- ScaleData(slide.seq)
		slide.seq <- RunPCA(slide.seq)
	}
	slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
	slide.seq <- FindClusters(slide.seq, resolution = resolution, verbose = F)
	assignments <- slide.seq@meta.data['seurat_clusters']
	colnames(assignments) <- 'cell_types'
	message(paste0(
		"cell types: ",
		paste(levels(assignments$cell_types), collapse = ', ')
	))
	if (use_silhouette) {
		distance_matrix <- dist(Embeddings(slide.seq[['pca']])[, 1:30])
	  clusters <- slide.seq@meta.data$seurat_clusters
	  silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
	  assignments$silhouette <- silhouette[, 3]

	  #assignments <- assignments[assignments$silhouette > 0, ]

	  cutoffs <- c()
	  for (cell_type in levels(assignments$cell_types)) {
	    cutoff <- quantile(
	    	assignments[assignments$cell_types == cell_type, ]$silhouette,
	    	probs=silhouette_cutoff
	    )
	    cutoffs <- append(cutoffs, cutoff)
	  }
	  names(cutoffs) <- levels(assignments$cell_types)
	  assignments <- assignments[
	  	assignments$silhouette > cutoffs[assignments$cell_types], ]
	}
	assignments
}


fit.gene.expression <- function(RCTD,
																cell_types = NULL,
																cell_type_threshold = 0,
																gene_threshold = 0) {
	if (is.null(cell_types))
		cell_types <- levels(RCTD@results$results_df$first_type)
	barcodes <- intersect(
		names(RCTD@spatialRNA@nUMI),
		colnames(RCTD@spatialRNA@counts)
	)
	X <- as.matrix(rep(1, length(barcodes)))
	rownames(X) <- barcodes
	RCTD <- run.CSIDE(
		RCTD, 
		X,
		barcodes,
		cell_types,
		cell_type_threshold = cell_type_threshold,
		gene_threshold = gene_threshold,
		sigma_gene = F,
		test_genes_sig = F,
		params_to_test = 1
	)
	cell_type_info <- list(
		as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)),
		cell_types,
		length(cell_types)
	)
	list(info = cell_type_info, renorm = cell_type_info)
}


assign.cell.types <- function(RCTD, assignments, weight = 0.9) {
	RCTD@internal_vars$cell_types_assigned <- T
	RCTD@config$RCTDmode <- "doublet"
	cell_type_names = levels(assignments[, 1])
	RCTD@cell_type_info <- list(
		info = list(data.frame(), cell_type_names, length(cell_type_names)),
		renorm = NULL
	)
	barcodes <- colnames(RCTD@spatialRNA@counts)
	N <- length(barcodes)
	weights_doublet = Matrix(0, nrow = N, ncol = 2)
	rownames(weights_doublet) = barcodes
	colnames(weights_doublet) = c('first_type', 'second_type')
	empty_cell_types = factor(character(N), levels = cell_type_names)
	spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_uncertain")
	results_df <- data.frame(
		spot_class = factor(character(N), levels=spot_levels),
		first_type = empty_cell_types,
		second_type = empty_cell_types
	)
	for(i in 1:N) {
		if(i %% 1000 == 0) print(paste("assign.cell.types: finished", i))
	  barcode <- barcodes[i]
		weights_doublet[i, ] = c(weight, 1 - weight)
		results_df[i, "spot_class"] = "singlet"
		results_df[i, "first_type"] = assignments[barcode, 1]
		results_df[i, "second_type"] = assignments[barcode, 1]
	}
	rownames(results_df) = barcodes
	RCTD@results <- list(
		results_df = results_df,
		weights_doublet = weights_doublet
	)
  RCTD
}


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
  RCTD
}