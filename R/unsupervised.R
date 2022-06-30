library(doParallel)
library(Matrix)
library(Seurat)
library(spacexr)


run.algorithm <- function(RCTD, 
                          doublet_mode = 'doublet',
                          cell_types = NULL,
                          n_iter = 50,
                          MIN_CHANGE = 0.001) {
  if (!(doublet_mode %in% c('doublet', 'full', 'subtype'))) {
    stop(paste0("run.unsupervised: doublet_mode=", doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, full, or subtype."))
  }
  RCTD@config$RCTDmode <- doublet_mode
  RCTD <- choose_sigma_c(RCTD)
  RCTD <- iter.optim(
    RCTD, 
    doublet_mode = doublet_mode, 
    cell_types = cell_types,
    n_iter = n_iter, 
    MIN_CHANGE = MIN_CHANGE
  )
}

iter.optim <- function(RCTD, 
                       doublet_mode = 'doublet',
                       cell_types = NULL,
                       n_iter = 50,
                       MIN_CHANGE = 1e-3) {
  if (is.null(cell_types)) {
    cell_types <- RCTD@cell_type_info$info[[2]]
  }
  barcodes <- intersect(
    names(RCTD@spatialRNA@nUMI), 
    colnames(RCTD@spatialRNA@counts)
  )
  X <- as.matrix(rep(1, length(barcodes)))
  rownames(X) <- barcodes
  RCTD@de_results$gene_fits$mean_val <- as.matrix(
    log(RCTD@cell_type_info$info[[1]])
  )
  originalSpatialRNA <- RCTD@originalSpatialRNA
  DOUBLET_THRESHOLD <- RCTD@config$DOUBLET_THRESHOLD
  RCTD@originalSpatialRNA <- RCTD@spatialRNA
  RCTD@config$DOUBLET_THRESHOLD <- 10
  RCTD@config$MIN_CHANGE_REG <- 1e-2
  RCTD@config$MIN_CHANGE_DE <- 1e-2

  message('iter.optim: assigning initial cell types')
  RCTD <- fitPixels(
    RCTD, 
    doublet_mode = doublet_mode
  )
  RCTD_list <- list(RCTD)

  for (i in 1:n_iter) {
    message(paste('iter.optim: running iteration', i))
    message('fitting gene expression profiles')
    initialSol <- RCTD@de_results$gene_fits$mean_val
    RCTD <- run.CSIDE(
      RCTD,
      X,
      barcodes,
      cell_types,
      doublet_mode = (doublet_mode == 'doublet'),
      cell_type_threshold = 0,
      gene_threshold = 0,
      sigma_gene = F,
      test_genes_sig = F,
      params_to_test = 1,
      logs = T,
      initialSol = initialSol[, cell_types]
    )
    info <- as.data.frame(exp(RCTD@de_results$gene_fits$mean_val))
    RCTD@cell_type_info$info[[1]][, cell_types] <- info 
    RCTD@cell_type_info$renorm[[1]][, cell_types] <- info 

    message('fitting cell types')
    if (doublet_mode == 'doublet') {
      initialSol <- RCTD@results$weights_pair
    } else {
      initialSol <- RCTD@results$weights
    }
    RCTD <- fitPixels(
      RCTD, 
      doublet_mode = doublet_mode,
      initialSol = initialSol
    )

    change <- weights_change(RCTD_list[[i]], RCTD)
    message(paste('change:', change))
    RCTD@config$MIN_CHANGE_REG <- max(min(1e-2, change ** 2), 1e-3)
    RCTD@config$MIN_CHANGE_DE <- max(min(1e-2, change ** 2), 1e-3)
    if (change < MIN_CHANGE) break
    RCTD_list[[i + 1]] <- RCTD
  }

  if (doublet_mode == 'doublet') {
    results <- RCTD@results$results_df
    reassign <- rownames(
      results[results$spot_class == 'doublet_certain' & 
      results$singlet_score - results$min_score < DOUBLET_THRESHOLD, ]
    )
    RCTD@results$results_df[reassign, ]$spot_class <- 'singlet'
  }
  RCTD@originalSpatialRNA <- originalSpatialRNA
  RCTD_list[[i + 1]] <- RCTD
  RCTD_list
}


initialize.clusters <- function(puck, resolution = 0.7, SCT = T) {
  message('Begin: initialize.clusters')
  slide.seq <- CreateSeuratObject(counts = puck@counts, assay = "Spatial")
  if (SCT) {
    slide.seq <- SCTransform(slide.seq, assay = "Spatial", verbose = F)
    slide.seq <- RunPCA(slide.seq, assay = "SCT", verbose = F)
  } else {
    slide.seq <- NormalizeData(slide.seq)
    slide.seq <- FindVariableFeatures(slide.seq)
    slide.seq <- ScaleData(slide.seq)
    slide.seq <- RunPCA(slide.seq)
  }
  slide.seq <- FindNeighbors(slide.seq, dims = 1:30, verbose = F)
  slide.seq <- FindClusters(slide.seq, resolution = resolution, verbose = F)
  clusters <- slide.seq@meta.data['seurat_clusters']
  colnames(clusters) <- 'cell_types'
  message(paste0("initialize.clusters: ", length(levels(clusters$cell_types)), " clusters generated"))
  cell_type_info_from_clusters(puck, clusters)
}

create.object <- function(spatialRNA,
                          cell_type_info,
                          max_cores = 4,
                          gene_cutoff_reg = 0.0002,
                          fc_cutoff_reg = 0.75,
                          UMI_min = 100,
                          UMI_max = 20000000, 
                          UMI_min_sigma = 300,
                          class_df = NULL,
                          gene_list_reg = NULL,
                          CONFIDENCE_THRESHOLD = 10,
                          DOUBLET_THRESHOLD = 25) {
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
    MIN_CHANGE_DE = 0.001,
    UMI_max = UMI_max,
    MIN_OBS = 3,
    CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD,
    DOUBLET_THRESHOLD = DOUBLET_THRESHOLD
  )
  reference <- new(
    "Reference",
    cell_types = factor(),
    counts = as(matrix(), 'dgCMatrix')
  )
  puck.original = restrict_counts(
    spatialRNA,
    rownames(spatialRNA@counts),
    UMI_thresh = config$UMI_min,
    UMI_max = config$UMI_max
  )
  if (is.null(gene_list_reg)) {
    message(paste("create.object: getting regression differentially expressed genes: "))
    gene_list_reg = get_de_genes(
      cell_type_info,
      puck.original,
      fc_thresh = config$fc_cutoff_reg,
      expr_thresh = config$gene_cutoff_reg,
      MIN_OBS = config$MIN_OBS
    )
  }
  if(length(gene_list_reg) == 0) {
    stop(paste("create.object: Error: 0 regression differentially expressed genes found"))
  }
  puck = restrict_counts(
    spatialRNA,
    gene_list_reg,
    UMI_thresh = config$UMI_min,
    UMI_max = config$UMI_max
  )
  puck = restrict_puck(puck, colnames(puck@counts))
  if (is.null(class_df)) {
    class_df <- data.frame(
      cell_type_info[[2]],
      row.names = cell_type_info[[2]]
    )
    colnames(class_df)[1] = "class"
  }
  cell_type_info <- list(info = cell_type_info, renorm = cell_type_info)
  cell_type_info$info[[1]] <- cell_type_info$info[[1]][gene_list_reg, ]
  cell_type_info$renorm[[1]] <- cell_type_info$renorm[[1]][gene_list_reg, ]
  internal_vars <- list(
    gene_list_reg = gene_list_reg,
    gene_list_bulk = gene_list_reg,
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

run.de.initialization <- function(RCTD, 
                                  doublet_mode = 'doublet',
                                  cell_types = NULL,
                                  n_iter = 50,
                                  MIN_CHANGE = 0.001) {
  cell_type_info <- cell_type_info_from_de(RCTD)
  RCTD <- create.object(
    RCTD@originalSpatialRNA,
    cell_type_info
  )
  RCTD <- run.algorithm(
    RCTD,
    doublet_mode = doublet_mode,
    cell_types = cell_types,
    n_iter = n_iter,
    MIN_CHANGE = MIN_CHANGE
  )
}

run.subtypes <- function(RCTD,
                         n_iter = 100,
                         MIN_CHANGE = 1e-3) {
  cell_types <- RCTD@internal_vars$subtypes
  RCTD@config$RCTDmode <- 'subtype'
  RCTD <- choose_sigma_c(RCTD)
  RCTD <- iter.optim(
    RCTD,
    doublet_mode = 'subtype',
    cell_types = cell_types,
    n_iter = n_iter,
    MIN_CHANGE = MIN_CHANGE
  )
}

initialize.subtypes <- function(RCTD,
                                cell_types, 
                                resolution = 0.7,
                                SCT = T,
                                gene_list = NULL,
                                fc_thresh = 0.5,
                                expr_thresh = 1e-4) {
  message('initialize.subtypes: gathering results')
  if (RCTD@config$RCTDmode == 'doublet') {
    weights <- weights_from_results(RCTD@results)
  } else if (RCTD@config$RCTDmode == 'full') {
    weights <- as.matrix(RCTD@results$weights)
  }
  weights <- weights[rowSums(as.matrix(weights[, cell_types])) > 0, ]
  cell_types_present <- apply(
    weights[, setdiff(colnames(weights), cell_types)], 
    1, 
    function(x) names(which(x > 0))
  )
  barcodes <- rownames(weights)
  singlet_barcodes <- barcodes[rowSums(as.matrix(weights[, cell_types])) == 1]
  singlet_puck <- restrict_puck(RCTD@originalSpatialRNA, singlet_barcodes)
  message('initialize.subtypes: getting subtype gene expression profiles: ')
  subtype_info <- initialize.clusters(
    singlet_puck,
    resolution = resolution,
    SCT = SCT
  )
  subtype_info[[2]] <- sapply(
    subtype_info[[2]], 
    function(x) paste0('subtype_', x)
  )
  colnames(subtype_info[[1]]) <- subtype_info[[2]]
  RCTD@spatialRNA <- RCTD@originalSpatialRNA
  if (is.null(gene_list)) {
    message('initialize.subtypes: getting subtype regression differentially expressed genes: ')
    subtype_genes <- get_de_genes(
      subtype_info, 
      singlet_puck, 
      fc_thresh = fc_thresh, 
      expr_thresh = expr_thresh, 
      MIN_OBS = 3
    )
    message('initialize.subtypes: getting supertype gene expression profiles: ')
    supertype_info <- cell_type_info_from_de(RCTD)
    message('initialize.subtypes: getting supertype highly expressed genes: ')
    supertype_genes <- get_de_genes(
      supertype_info, 
      RCTD@originalSpatialRNA, 
      cell_types = cell_types,
      fc_thresh = 0, 
      expr_thresh = 1e-4, 
      MIN_OBS = 3
    )
    gene_list_tot <- union(subtype_genes, supertype_genes)
  } else {
    message('initialize.subtypes: getting supertype gene expression profiles: ')
    supertype_info <- cell_type_info_from_de(RCTD, gene_list = gene_list)
    gene_list_tot <- gene_list
  }
  info <- cbind(
    supertype_info[[1]][gene_list_tot, ], 
    subtype_info[[1]][gene_list_tot, ]
  )
  info <- info[, setdiff(colnames(info), cell_types)]
  cell_type_info <- list(info, colnames(info), length(colnames(info)))
  RCTD = create.object(
    restrict_puck(RCTD@originalSpatialRNA, barcodes), 
    cell_type_info, 
    gene_list_reg = gene_list_tot
  )
  cell_types_present <- lapply(
    cell_types_present, 
    function(x) c(x, subtype_info[[2]])
  )
  RCTD@internal_vars$subtypes = subtype_info[[2]]
  RCTD@internal_vars$cell_types_present = cell_types_present
  RCTD
}


weights_from_results <- function(results) {
  df <- results$results_df
  wd <- results$weights_doublet
  singlets <- rownames(df[df$spot_class == 'singlet', ])
  doublets <- rownames(df[df$spot_class == 'doublet_certain' |
                          df$spot_class == 'doublet_uncertain', ])
  barcodes <- c(singlets, doublets)
  cell_types <- levels(df$first_type)
  weights <- matrix(0, nrow = length(barcodes), ncol = length(cell_types))
  rownames(weights) <- barcodes; colnames(weights) <- cell_types
  for (b in singlets) {
    weights[b, as.character(df[b, 'first_type'])] <- 1
  }
  for (b in doublets) {
    weights[b, as.character(df[b, 'first_type'])] <- wd[b, 'first_type']
    weights[b, as.character(df[b, 'second_type'])] <- wd[b, 'second_type']
  }
  weights
}

cell_type_info_from_clusters <- function(puck, clusters) {
  counts <- puck@counts[, rownames(clusters)]
  nUMI <- puck@nUMI[rownames(clusters)]
  cell_types <- clusters$cell_types
  names(cell_types) <- rownames(clusters)
  get_cell_type_info(counts, cell_types, nUMI)
}

cell_type_info_from_singlets <- function(RCTD,
                                         cell_types = NULL) {
  results <- RCTD@results$results_df
  if (is.null(cell_types)) {
    cell_types <- levels(results$first_type)
  }
  singlets <- results[results$spot_class == 'singlet',]
  cell_type_list <- singlets$first_type
  names(cell_type_list) <- rownames(singlets)
  get_cell_type_info(
    RCTD@originalSpatialRNA@counts[, rownames(singlets)],
    cell_type_list,
    RCTD@originalSpatialRNA@nUMI[rownames(singlets)],
    cell_type_names = cell_types
  )
}

cell_type_info_from_de <- function(RCTD,
                                   cell_types = NULL,
                                   gene_list = NULL,
                                   cell_type_threshold = 0,
                                   gene_threshold = 0) {
  if (is.null(cell_types)) {
    cell_types <- levels(RCTD@results$results_df$first_type)
  }
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
    gene_list_tot = gene_list,
    doublet_mode = (RCTD@config$RCTDmode == "doublet"),
    cell_type_threshold = cell_type_threshold,
    gene_threshold = gene_threshold,
    sigma_gene = F,
    test_genes_sig = F,
    params_to_test = 1,
    logs = T
  )
  list(
    as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)),
    cell_types,
    length(cell_types)
  )
}

process_cell_type_info <- function(reference, cell_type_names, CELL_MIN = 25) {
   message("Begin: process_cell_type_info")
   message(paste("process_cell_type_info: number of cells in reference:", dim(reference@counts)[2]))
   message(paste("process_cell_type_info: number of genes in reference:", dim(reference@counts)[1]))
   cell_counts = table(reference@cell_types)
   print(cell_counts)

   if(min(cell_counts) < CELL_MIN)
      stop(paste0("process_cell_type_info error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
   cell_type_info <- get_cell_type_info(reference@counts, reference@cell_types, reference@nUMI
                                        , cell_type_names = cell_type_names)
   message("End: process_cell_type_info")
   return(cell_type_info)
}


weights_change <- function(RCTD1, RCTD2) {
  weights1 <- RCTD1@results$weights
  weights2 <- RCTD2@results$weights
  norm(weights1 - weights2) / dim(weights1)[1]
}

get_marker_data <- function(cell_type_names, cell_type_means, gene_list) {
  marker_means = cell_type_means[gene_list,]
  marker_norm = marker_means / rowSums(marker_means)
  marker_data = as.data.frame(cell_type_names[max.col(marker_means)])
  marker_data$max_epr <- apply(cell_type_means[gene_list,],1,max)
  colnames(marker_data) = c("cell_type",'max_epr')
  rownames(marker_data) = gene_list
  marker_data$log_fc <- 0
  epsilon <- 1e-9
  for(cell_type in unique(marker_data$cell_type)) {
    cur_genes <- gene_list[marker_data$cell_type == cell_type]
    other_mean = rowMeans(cell_type_means[cur_genes,cell_type_names != cell_type])
    marker_data$log_fc[marker_data$cell_type == cell_type] <- log(epsilon + cell_type_means[cur_genes,cell_type]) - log(epsilon + other_mean)
  }
  return(marker_data)
}

plot_all_cell_types <- function(results_df, coords, cell_type_names, resultsdir, size=0.15) {
  barcodes = rownames(results_df[results_df$spot_class != "reject" & results_df$first_type %in% cell_type_names,])
  my_table = coords[barcodes,]
  my_table$class = results_df[barcodes,]$first_type
  n_levels = length(levels(my_table$class))
  my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
  pres = unique(as.integer(my_table$class))
  pres = pres[order(pres)]
  if(n_levels > 21)
    my_pal = pals::polychrome(n_levels)
  if(n_levels > 36)
    stop("Plotting currently supports at most 36 cell types as colors")
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = size, shape=19,color=class)) +
    ggplot2::scale_color_manual(values = my_pal[pres])+ ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
  pdf(file.path(resultsdir,"all_cell_types.pdf"))
  invisible(print(plot))
  dev.off()
  return(plot)
}

singlet_table <- function(truth_results, pred_results) {
  truth_singlet <- truth_results[truth_results$spot_class == 'singlet',]
  pred_singlet <- pred_results[pred_results$spot_class == 'singlet',]
  common_barcode <- intersect(row.names(truth_singlet), row.names(pred_singlet))
  truth_singlet <- truth_singlet[common_barcode,]
  pred_singlet <- pred_singlet[common_barcode,]
  truth_types = unlist(list(truth_singlet[,'first_type']))
  pred_types = unlist(list(pred_singlet[,'first_type']))
  return(table(truth_types, pred_types))
}


# not used anymore

weights_mae <- function(RCTD1, RCTD2) {
  weights1 <- RCTD1@results$weights
  weights2 <- RCTD2@results$weights
  sum(abs(weights1 - weights2)) / length(weights1)
}


assignment_accuracy <- function(RCTD1, RCTD2) {
  results1 <- RCTD1@results$results_df[, 1:3]
  results2 <- RCTD2@results$results_df[, 1:3]
  placeholder <- levels(results1$first_type)[1]
  if (length(results1[results1$spot_class == 'singlet', ]$second_type) > 0)
    results1[results1$spot_class == 'singlet', ]$second_type <- placeholder
  if (length(results1[results1$spot_class == 'reject', ]$first_type) > 0)
    results1[results1$spot_class == 'reject', ]$first_type <- placeholder
  if (length(results1[results1$spot_class == 'reject', ]$second_type) > 0)
    results1[results1$spot_class == 'reject', ]$second_type <- placeholder
  if (length(results2[results2$spot_class == 'singlet', ]$second_type) > 0)
    results2[results2$spot_class == 'singlet', ]$second_type <- placeholder
  if (length(results2[results2$spot_class == 'reject', ]$first_type) > 0)
    results2[results2$spot_class == 'reject', ]$first_type <- placeholder
  if (length(results2[results2$spot_class == 'reject', ]$second_type) > 0)
    results2[results2$spot_class == 'reject', ]$second_type <- placeholder
  common <- as.data.frame(results1 == results2)
  correct <- dim(common[common$spot_class & common$first_type & common$second_type, ])[1]
  correct / dim(results1)[1]
}

normalize_reference <- function(puck, reference, cell_type_names = NULL) {
  if (is.null(cell_type_names)) {
    cell_type_names <- levels(reference@cell_types)
  }
  cell_type_info = list(
    info = process_cell_type_info(
      reference,
      cell_type_names = cell_type_names,
      CELL_MIN = 25
    ),
    renorm = NULL
  )
  puck = restrict_counts(
    puck,
    rownames(puck@counts),
    UMI_thresh = 100,
    UMI_max = 20000000
  )
  gene_list_bulk = intersect(rownames(puck@counts), rownames(reference@counts))
  RCTD <- new(
    "RCTD",
    spatialRNA = puck,
    originalSpatialRNA = puck,
    reference = reference,
    config = list(MIN_CHANGE_BULK = 0.0001), 
    cell_type_info = cell_type_info, 
    internal_vars = list(gene_list_bulk = gene_list_bulk, proportions = NULL)
  )
  RCTD <- fitBulk(RCTD)
  RCTD@cell_type_info
}

# semisupervised learning, does not work right now

initialize.semisupervised <- function(puck,
                                      reference, 
                                      cell_types = NULL, 
                                      resolution = 0.7, 
                                      SCT = T) {
  message(paste0("initialize.semisupervised: generating clusters"))
  clusters <- initialize.clusters(
    puck, 
    resolution = resolution, 
    SCT = SCT
  )
  message(paste0("initialize.semisupervised: normalizing reference"))
  normref <- normalize_reference(
    puck, 
    reference, 
    cell_type_names = cell_types
  )
  normref$renorm[[1]][is.na(normref$renorm[[1]])] <- 0
  clusters$renorm[[1]] <- clusters$renorm[[1]][rownames(normref$renorm[[1]]), ]
  gene_cor <- cor(clusters$renorm[[1]], normref$renorm[[1]])
  clusters$renorm[[1]] <- clusters$renorm[[1]][, -apply(gene_cor, 2, which.max)]
  message(paste0("initialize.semisupervised: ", dim(clusters$renorm[[1]])[2], " learned cell type(s)"))
  info <- cbind(normref$renorm[[1]], clusters$renorm[[1]])
  info <- list(info, colnames(info), length(colnames(info)))
  list(info = info, renorm = info)
}


initialize.semisupervised.2 <- function(puck,
                                        reference, 
                                        cell_types = NULL, 
                                        resolution = 0.7, 
                                        SCT = T) {
  message(paste0("initialize.semisupervised: generating clusters"))
  clusters <- initialize.clusters(
    puck, 
    resolution = resolution, 
    SCT = SCT
  )
  message(paste0("initialize.semisupervised: processing reference"))
  normref <- list(
    info = process_cell_type_info(
      reference,
      cell_type_names = cell_types,
      CELL_MIN = 25
    ),
    renorm = NULL
  )
  normref$renorm <- normref$info
  gene_list <- intersect(rownames(puck@counts), rownames(reference@counts))
  clusters$renorm[[1]] <- clusters$renorm[[1]][gene_list, ]
  normref$renorm[[1]] <- normref$renorm[[1]][gene_list, ]
  gene_cor <- cor(clusters$renorm[[1]], normref$renorm[[1]])
  clusters$renorm[[1]] <- clusters$renorm[[1]][, -apply(gene_cor, 2, which.max)]
  message(paste0("initialize.semisupervised: ", dim(clusters$renorm[[1]])[2], " learned cell type(s)"))
  info <- cbind(normref$renorm[[1]], clusters$renorm[[1]])
  info <- list(info, colnames(info), length(colnames(info)))
  list(info = info, renorm = info)
}