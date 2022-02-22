library(Seurat)
library(cluster)


#' Generates Seurat cell type clusters of a SpatialRNA object.
#' 
#' @param puck a SpatialRNA object.
#' @param resolution (default 0.3) the resolution of the clustering. A higher resolution results in more clusters.
#' @param SCT (default TRUE) whether to use regularized negative binomial regression to normalize count data.
#' @return a dataframe with bead barcodes as rownames and assigned cluster number as values.
#' @export
gen.clusters <- function(puck, resolution = 1, SCT = T, silhouette_cutoff = 0) {
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
  message(paste0("cell types: ",paste(levels(assignments$cell_types), collapse = ', ')))
  # assign silhouette scores
  distance_matrix <- dist(Embeddings(slide.seq[['pca']])[, 1:30])
  clusters <- slide.seq@meta.data$seurat_clusters
  silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
  assignments$silhouette <- silhouette[, 3]
  assignments <- assignments[assignments$silhouette > silhouette_cutoff,]
  return(assignments)
}


# plot silhouette score
slide.seq@meta.data$silhouette_score <- silhouette[,3]
mean_silhouette_score <- mean(slide.seq@meta.data$silhouette_score)
p <- slide.seq@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  arrange(seurat_clusters,-silhouette_score) %>%
  mutate(barcode = factor(barcode, levels = barcode)) %>%
  ggplot() +
  geom_col(aes(barcode, silhouette_score, fill = seurat_clusters), show.legend = FALSE) +
  geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
  scale_x_discrete(name = 'Cells') +
  scale_y_continuous(name = 'Silhouette score') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
