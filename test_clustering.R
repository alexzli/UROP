library(Seurat)

puck <- readRDS('../UROP/data/puckCropped_cerebellum_slideseq.rds')
slide.seq <- CreateSeuratObject(counts = puck@counts)

slide.seq <- NormalizeData(slide.seq)
slide.seq <- FindVariableFeatures(slide.seq)
slide.seq <- ScaleData(slide.seq)
#slide.seq <- SCTransform(slide.seq, ncells = 3000, verbose = FALSE)

slide.seq <- RunPCA(slide.seq)
#slide.seq <- RunTSNE(slide.seq, dims = 1:30)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq.clusters <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
plot1 <- DimPlot(slide.seq.clusters, reduction = "umap", label = TRUE)
plot1
