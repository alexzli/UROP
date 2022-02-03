library(caret)
library(doParallel)
library(Matrix)
library(Seurat)
library(spacexr)

run.semisupervised <- function(RCTD, max_cores = 4, keep_gene_list = T) {
	if (keep_gene_list) {
		gene_list <- RCTD@internal_vars$gene_list_reg
		RCTD@originalSpatialRNA <- restrict_counts(RCTD@originalSpatialRNA, gene_list, UMI_thresh = 100, UMI_max = 20000000)
		gene_threshold = 0
	} else {
		gene_list <- NULL
		gene_threshold = 5e-5
	}
	cell_type_info <- fit.gene.expression(RCTD, gene_threshold = gene_threshold)
	myRCTD <- create.RCTD.noref(RCTD@originalSpatialRNA, max_cores = max_cores, cell_type_info = cell_type_info, 
															gene_list_reg = gene_list, class_df = RCTD@internal_vars$class_df)
	return(iter.optim(myRCTD))
}

run.unsupervised <- function(puck, max_cores = 4, resolution = 1, SCT = F, info_type = 'mean') {
	assignments <- gen.clusters(puck, resolution = resolution, SCT = SCT)
	if (info_type == 'mean') {
		cell_type_info <- cell_type_info_from_assignments(puck, assignments)
	} else if (info_type == 'de') {
		cell_type_info <- de_info_from_assignments(puck, assignments, gene_threshold = 5e-5)
	} else {
		stop(paste0("run.unsupervised: info_type=",doublet_mode, " is not a valid choice. Please set info_type=mean or de."))
	}
	myRCTD <- create.RCTD.noref(puck, max_cores = max_cores, cell_type_info = cell_type_info)
	return(iter.optim(myRCTD))
}

gen.clusters <- function(puck, resolution = 1, SCT = F) {
	slide.seq <- CreateSeuratObject(counts = puck@counts)
	# data preprocessing
	if (SCT) {
		slide.seq <- SCTransform(slide.seq, ncells = 3000, verbose = F)
	} else {
		slide.seq <- NormalizeData(slide.seq)
		slide.seq <- FindVariableFeatures(slide.seq)
		slide.seq <- ScaleData(slide.seq)
	}
	# dimensionality reduction and clustering
	slide.seq <- RunPCA(slide.seq)
	slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
	slide.seq <- FindClusters(slide.seq, resolution = resolution, verbose = F)
	assignments <- slide.seq@meta.data['seurat_clusters']
	colnames(assignments) <- 'cell_types'
	return(assignments)
}

cell_type_info_from_assignments <- function(puck, assignments) {
	assigned_cell_types <- assignments[colnames(puck@counts),]
	names(assigned_cell_types) <- rownames(assignments)
	info <- get_cell_type_info(puck@counts, assigned_cell_types, puck@nUMI)
	return(list(info = info, renorm = info))
}

create.RCTD.noref <- function(spatialRNA, max_cores = 4, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 20000000, 
							  UMI_min_sigma = 300, MAX_MULTI_TYPES = 4, cell_type_info = NULL, gene_list_reg = NULL, class_df = NULL) {

	config <- list(gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma, max_cores = max_cores,
	               N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30, MIN_CHANGE_BULK = 0.0001, MIN_CHANGE_REG = 0.001, UMI_max = UMI_max, MIN_OBS = 3, MAX_MULTI_TYPES = MAX_MULTI_TYPES)
	reference <- new("Reference", cell_types = factor(), counts = as(matrix(), 'dgCMatrix'))
	puck.original = restrict_counts(spatialRNA, rownames(spatialRNA@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
	if (is.null(cell_type_info)) {
		cell_type_info <- list(info = list(data.frame(), c(), 0), renorm = list(data.frame(), c(), 0))
		puck = puck.original
	} else {
		message("create.RCTD.noref: getting regression differentially expressed genes: ")
		if (is.null(gene_list_reg))
			gene_list_reg = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
		if(length(gene_list_reg) == 0)
		  stop("create.RCTD.noref: Error: 0 regression differentially expressed genes found")
		puck = restrict_counts(spatialRNA, gene_list_reg, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
		puck = restrict_puck(puck, colnames(puck@counts))
		puck.original = puck
		if (is.null(class_df))
			class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]]); colnames(class_df)[1] = "class"
		cell_type_info$info[[1]] <- cell_type_info$info[[1]][gene_list_reg,]; cell_type_info$renorm[[1]] <- cell_type_info$renorm[[1]][gene_list_reg,]
	}
	internal_vars <- list(gene_list_reg = gene_list_reg, proportions = NULL, class_df = class_df, cell_types_assigned = F)
	new("RCTD", spatialRNA = puck, originalSpatialRNA = puck.original, reference = reference, config = config, cell_type_info = cell_type_info, internal_vars = internal_vars)
}

iter.optim <- function(RCTD, cell_types = NULL, CELL_MIN_INSTANCE = 0, max_n_iter = 20, convergence_cutoff = 0.999, discovery_threshold = 0.99) {
	RCTD@config$RCTDmode <- 'doublet'
	RCTD <- choose_sigma_c(RCTD)
	message('run.iter.optim: assigning initial cell types')
	RCTD <- fitPixels(RCTD, doublet_mode = 'doublet')
	barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
	if (is.null(cell_types)) {
		cell_type_count <- aggregate_cell_types(RCTD, barcodes, doublet_mode = T)
		cell_types <- names(which(cell_type_count >= CELL_MIN_INSTANCE))
	}
	RCTD_list = list()
	RCTD_list[[1]] <- RCTD
	X <- as.matrix(rep(1, length(barcodes)))
	rownames(X) <- barcodes
	for (i in 1:max_n_iter) {
		message(paste('run.iter.optim: running iteration', i))
		message('fitting gene expression profiles')
		RCTD <- run.CSIDE(RCTD, X, barcodes, cell_types, cell_type_threshold = 0, gene_threshold = 0, sigma_gene = F, test_genes_sig = F, params_to_test = 1)
		message('fitting cell types')
		info <- list(as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)), cell_types, length(cell_types))
		cell_type_info <- list(info = info, renorm = info)
		RCTD@cell_type_info <- cell_type_info
		RCTD <- fitPixels(RCTD)
		RCTD_list[[i+1]] <- RCTD
		accuracy <- singlet_accuracy(RCTD_list[[i]], RCTD_list[[i+1]], cell_types)
		n_singlets <- dim(RCTD@results$results_df[RCTD@results$results_df$spot_class == 'singlet',])[1]
		prev_singlets <- dim(RCTD_list[[i]]@results$results_df[RCTD_list[[i]]@results$results_df$spot_class == 'singlet',])[1]
		message(paste('classification convergence accuracy:', accuracy))
		if (accuracy >= convergence_cutoff & discovery_threshold <= n_singlets/prev_singlets & 1/discovery_threshold >= n_singlets/prev_singlets)
			break
	}
	return(RCTD_list)
}

singlet_accuracy <- function(RCTD_ref, RCTD_pred, cell_types) {
	ref_singlet <- RCTD_ref@results$results_df
	ref_singlet <- ref_singlet[ref_singlet$spot_class == 'singlet' & is.element(ref_singlet$first_type, cell_types),]
	pred_singlet <- RCTD_pred@results$results_df
	pred_singlet <- pred_singlet[pred_singlet$spot_class == 'singlet' & is.element(pred_singlet$first_type, cell_types),]
	common_barcode <- intersect(row.names(ref_singlet), row.names(pred_singlet))
	ref_singlet <- ref_singlet[common_barcode,]
	pred_singlet <- pred_singlet[common_barcode,]
	ref_types = unlist(list(ref_singlet[,'first_type']))
	pred_types = unlist(list(pred_singlet[,'first_type']))
	return(confusionMatrix(pred_types, ref_types)$overall['Accuracy'])
}

de_info_from_assignments <- function(puck, assignments, gene_threshold = 0) {
	myRCTD <- create.RCTD.noref(puck)
	myRCTD <- set_internal_vars(myRCTD)
	myRCTD <- assign.cell.types(myRCTD, assignments)
	return(fit.gene.expression(myRCTD, gene_threshold = gene_threshold))
}

assign.cell.types <- function(RCTD, assignments, weight = 0.9) {
	RCTD@internal_vars$cell_types_assigned <- T
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

fit.gene.expression <- function(RCTD, cell_types = NULL, cell_type_threshold = 0, gene_threshold = 0) {
	if (is.null(cell_types))
		cell_types <- levels(RCTD@results$results_df$first_type)
	barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
	X <- rep(1, length(barcodes))
	X <- as.matrix(X)
	rownames(X) <- barcodes
	RCTD <- run.CSIDE(RCTD, X, barcodes, cell_types, cell_type_threshold = cell_type_threshold, gene_threshold = gene_threshold, 
										sigma_gene = F, test_genes_sig = F, params_to_test = 1)
	cell_type_info <- list(as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)), cell_types, length(cell_types))
	return(list(info = cell_type_info, renorm = cell_type_info))
}