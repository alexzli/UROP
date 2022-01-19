library(Seurat)

#' Generates Seurat cell type clusters of a SpatialRNA object.
#' 
#' @param puck a SpatialRNA object.
#' @param resolution (default 0.3) the resolution of the clustering. A higher resolution results in more clusters.
#' @param SCT (default FALSE) whether to use regularized negative binomial regression to normalize count data.
#' @return a dataframe with bead barcodes as rownames and assigned cluster number as values.
#' @export
gen.clusters <- function(puck, resolution = 0.3, SCT = F) {
	slide.seq <- CreateSeuratObject(counts = puck@counts)
	# data preprocessing
	if (SCT)
		slide.seq <- SCTransform(slide.seq, ncells = 3000, verbose = FALSE)
	else {
		slide.seq <- NormalizeData(slide.seq)
		slide.seq <- FindVariableFeatures(slide.seq)
		slide.seq <- ScaleData(slide.seq)
	}
	# dimensionality reduction and clustering
	slide.seq <- RunPCA(slide.seq)
	slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
	slide.seq <- FindClusters(slide.seq, resolution = resolution, verbose = FALSE)
	assignments <- slide.seq@meta.data['seurat_clusters']
	colnames(assignments) <- 'cell_types'
	return(assignments)
}
