library(Seurat)

puck <- readRDS('../puckCropped.rds')
slide.seq <- CreateSeuratObject(counts = puck@counts)

slide.seq <- NormalizeData(slide.seq)
slide.seq <- FindVariableFeatures(slide.seq)
slide.seq <- ScaleData(slide.seq)
#slide.seq <- SCTransform(slide.seq, ncells = 3000, verbose = FALSE)

slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot1