
# Code for generating MERFISH dataset as used in the STdeconovlve paper.


## Moffit et al. 2018 raw data downloaded from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248/
moffit <- read.csv2(file = "../UROP/data/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv",
                    sep = ",")
## split into metadata
annot.table_ <- moffit[,c(1:9)]
## and gene counts
counts_ <- moffit[,-c(1:9)]

## convert gene counts to numerics
counts_ <- as.matrix(data.frame(apply(counts_, 2, function(x) as.numeric(as.character(x)))))

rownames(annot.table_) <- annot.table_[,"Cell_ID"]
rownames(counts_) <- annot.table_[,"Cell_ID"]

## Moffit et al. 2018 Supplemental Table S6 contains the genes profiled
merfish_genes <- readxl::read_excel(path = "../UROP/data/aau5324_Moffitt_Table-S6.xlsx",
                                    range = "A2:D163")
## find the smFISH genes
smFISH_genes <- dplyr::filter(merfish_genes, Barcode == "Sequential stain")

## filter for MERFISH profiled genes
counts_ <- counts_[, which(!colnames(counts_) %in% smFISH_genes$`Gene name`)]
## remove Blank measurements
counts_ <- counts_[, !grepl("Blank", colnames(counts_))]

## transform the expression values to discrete counts
## divide by 1000 and multiply by each cell's volume

## "171023_FN7_1_M22_M26" (posterior) and "171021_FN7_2_M22_M26" (anterior) are datasets for Female Naive "Animal_ID" = 2.
## Together they account for all 12 bregma tissue sections.

## collect the pertinent metadata for the cells in the tissue sections
annot.table_ <- annot.table_[annot.table_$Animal_ID == 2,
                             c('Centroid_X', 'Centroid_Y', 'Bregma', "Cell_class", "Neuron_cluster_ID")]
#annot.table_ <- annot.table_[annot.table_$Animal_ID == 1 & annot.table_$Bregma == 0.01,
#                             c('Centroid_X', 'Centroid_Y', 'Bregma', "Cell_class", "Neuron_cluster_ID")]

## collapse the Oligodendrocytes and Endothelial cell-types
annot.table_[grep(pattern = "OD Mature",
                                x = annot.table_$Cell_class),]$Cell_class <- "OD Mature"
annot.table_[grep(pattern = "OD Immature",
                                x = annot.table_$Cell_class),]$Cell_class <- "OD Immature"
annot.table_[grep(pattern = "Endothelial",
                                x = annot.table_$Cell_class),]$Cell_class <- "Endothelial"

## remove "Ambiguous" cell-types
annot.table_ <- annot.table_[annot.table_$Cell_class != "Ambiguous",]
dim(annot.table_)

## remove rows with NA
annot.table_ <- na.omit(annot.table_)

## change the type of some of the metadata columns to numerics
annot.table_$Centroid_X <- as.numeric(as.character(annot.table_$Centroid_X))
annot.table_$Centroid_Y <- as.numeric(as.character(annot.table_$Centroid_Y))

## filter the gene counts matrix to only include the filtered cells
counts_ <- counts_[rownames(annot.table_),]
dim(counts_)

## retrieve raw counts
counts_ <- counts_/apply(counts_, 1, function(x) min(x[x > 0]))
counts_ <- apply(counts_, c(1,2), as.integer)

library(STdeconvolve)

## get list of major cell-types of each cell
majorCellTypes <- annot.table_$Cell_class
length(majorCellTypes)

## get a similar list but with neuronal cell-types expanded to subtypes
neuronalCellsubtypes <- unlist(lapply(rownames(annot.table_), function(cell){
  class <- annot.table_[cell,]$Cell_class
  neuron <- annot.table_[cell,]$Neuron_cluster_ID
  if (neuron != ""){
    i <- neuron
  } else {
    i <- class
  }
  i
}))
length(neuronalCellsubtypes)

## Using the gene counts of the individual cells, generate a ground truth gene expression profile
## for each of the cell-types

## major cell-types:
cellTypes <- majorCellTypes
cells <- rownames(annot.table_)
mat <- counts_[cells,]
mm <- model.matrix(~ 0 + factor(cellTypes))
colnames(mm) <- levels(factor(cellTypes))
gtGexpCellTypes <- t(t(as.matrix(mat)) %*% mm)
gtGexpCellTypes <- gtGexpCellTypes/rowSums(gtGexpCellTypes)
dim(gtGexpCellTypes)

## expanded neuronal cell-types:
cellTypes <- neuronalCellsubtypes
cells <- rownames(annot.table_)
mat <- counts_[cells,]
mm <- model.matrix(~ 0 + factor(cellTypes))
colnames(mm) <- levels(factor(cellTypes))
gtGexpNeuronalCellTypes <- t(t(as.matrix(mat)) %*% mm)
gtGexpNeuronalCellTypes <- gtGexpNeuronalCellTypes/rowSums(gtGexpNeuronalCellTypes)
dim(gtGexpNeuronalCellTypes)

## The function takes into account the cell-types in the "Cell_class" column of the metadata
## so make sure this is set to the cell-types of interest
annot.table_$Cell_class <- majorCellTypes

# updated to include nUMI in cellTypeTable
simulateBregmaSpots <- function (cellCentroidsAndClass, counts, patch_size = 100) {
  
  # dictionary hash table
  h <- hash::hash()
  
  data <- cellCentroidsAndClass
  
  for (bregma in unique(data$Bregma)) {
    
    bregma_key <- as.character(bregma)
    print(bregma_key)
    
    selected_bregma <- data[which(data$Bregma == bregma),]
    
    # 1. Get patch edge coordinates:
    
    # Sequence of X-coord positions for left edge of each patch:
    x_edges <- seq(min(selected_bregma$Centroid_X), max(selected_bregma$Centroid_X), patch_size)
    # drop first and last to avoid any issues with the edges of the whole region
    inner_x_edges <- x_edges[2:(length(x_edges)-1)]
    # Sequence of Y-coord positions for bottom edge of each patch:
    y_edges <- seq(min(selected_bregma$Centroid_Y), max(selected_bregma$Centroid_Y), patch_size)
    inner_y_edges <- y_edges[2:(length(y_edges)-1)]
    
    selected_bregma$patch_id <- character(length(rownames(selected_bregma)))
    
    # 2. add patch IDs to bregma cells, for the patch they belong to:
    
    for (x in inner_x_edges) {
      for (y in inner_y_edges) {
        patch_id <- paste0(as.character(x), "_", as.character(y))
        patch <- selected_bregma[which( (selected_bregma$Centroid_X > x) &
                                          (selected_bregma$Centroid_X < x+patch_size) &
                                          (selected_bregma$Centroid_Y > y) &
                                          (selected_bregma$Centroid_Y < y+patch_size) ),]
        
        if (length(rownames(patch)) > 0) {
          selected_bregma[rownames(patch),]$patch_id <- patch_id
        }
      }
    }
    
    # 3. get table of counts of cell types for each patch in bregma
    selected_bregma$nUMI <- rowSums(counts[rownames(selected_bregma),])
    patch_cell_types <- selected_bregma[which(selected_bregma$patch_id != ""),
                                        c("patch_id", "Cell_class", "nUMI")]
    patch_cell_types <- as.data.frame(lapply(patch_cell_types, rep, patch_cell_types$nUMI))
    patch_cell_types <- patch_cell_types[,c("patch_id","Cell_class")]
    cellTypeTable <- table(patch_cell_types[])
    
    # 4. total cell counts in each patch
    patchTotalCells <- rowSums(cellTypeTable)
    
    # 5. counts of unique cell types in each patch
    cellTypeCount <- c()
    for (i in seq_len(length(rownames(cellTypeTable)))) {
      patch_num_cell_types <- length(which(cellTypeTable[i,] != 0))
      cellTypeCount <- append(cellTypeCount, patch_num_cell_types)
    }
    
    # 6. collapse gene counts for cells in same spot to make simulation
    patches <- unique(selected_bregma$patch_id[which(!selected_bregma$patch_id == "")])
    patchGexp <- do.call(rbind, lapply(patches, function(patch){
      cells <- rownames(selected_bregma[which(selected_bregma$patch_id == patch),])
      mat <- as.matrix(counts[cells,])
      if (length(cells) == 1){
        patch_counts <- as.vector(mat)
      } else if (length(cells) > 1){
        patch_counts <- colSums(mat)
      } else if (length(cells) == 0){
        cat("WARNING:", bregma, "patch", patch, "had no cells in `counts` to use for simulated gene expression", "\n")
      }
      patch_counts
    }))
    rownames(patchGexp) <- patches
    
    # 7. gene count matrix for individual cells in the bregma
    # this also includes cells not assigned to patches
    bregma_cells <- rownames(selected_bregma)
    cellGexp <- counts[bregma_cells,]
    
    # 8. combine data objects and append to hash table
    h[[bregma_key]] <- list(bregmaFullDf = selected_bregma,
                            cellTypeTable = cellTypeTable,
                            patchTotalCells = patchTotalCells,
                            cellTypeCount = cellTypeCount,
                            cellGexp = cellGexp,
                            patchGexp = patchGexp)
  }
  return(h)
}

FN7_hash <- simulateBregmaSpots(annot.table_,
                                counts = counts_, # gene counts matrix
                                patch_size = 100) # size of the simulated pixels in um2 (units of the Centroids)

# table of number of cell types in each spot for entire FN7
# spot IDs for all spots in FN7
FN7_spotIDs <- unlist(lapply(hash::keys(FN7_hash), function(ix){
  # table to df
  rownames(FN7_hash[[ix]]$cellTypeTable)
}))
# combine all the cellTypeTables of each bregma in FN7 (counts of each cell type in each spot)
FN7_cellTypeTable <- lapply(hash::keys(FN7_hash), function(ix){
  # table to df
  as.data.frame.matrix(FN7_hash[[ix]]$cellTypeTable)
})
# combine into single df, and because some bregmas may be missing cell types,
# use rbindlist to keep all columns and add NAs to spots for cell types
# they are missing
FN7_cellTypeTable <- data.table::rbindlist(FN7_cellTypeTable, fill = TRUE)
# replace NAs with 0s
FN7_cellTypeTable[is.na(FN7_cellTypeTable)] <- 0
# spot IDs as row names
FN7_cellTypeTable <- as.matrix(FN7_cellTypeTable)
rownames(FN7_cellTypeTable) <- FN7_spotIDs

simBregmasFN7 <- lapply(hash::keys(FN7_hash), function(ix){
  bregma <- buildBregmaCorpus(hashTable = FN7_hash, 
                                bregmaID = ix)
  print(bregma$sim)
  print(dim(bregma$gtSpotTopics))
  print(dim(bregma$gtCtGenes))
  bregma
})

names(simBregmasFN7) <- hash::keys(FN7_hash)

# 1. sim
# combine the sim slam matrices to make the corpus for all spots across all bregmas
sim_N7 <- slam::as.simple_triplet_matrix(do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- as.matrix(simBregmasFN7[[ix]]$sim)
  m
})))
sim_N7

# 2. gtSpotTopics
# each gtSpotTopics ref can have different numbers of cell types that are present in each bregma,
# at least for the neuro. So need a way to combine and have all 75 columns,
# and set 0 for patches that do not have any of one ct
gtSpotTopics_N7 <- lapply(names(simBregmasFN7), function(ix){
  simBregmasFN7[[ix]]$gtSpotTopics
})
names(gtSpotTopics_N7) <- names(simBregmasFN7)
gtSpotTopics_N7 <- data.table::rbindlist(gtSpotTopics_N7, fill = TRUE)
dim(gtSpotTopics_N7)

gtSpotTopics_N7 <- as.matrix(gtSpotTopics_N7)
rownames(gtSpotTopics_N7) <- rownames(sim_N7)
# Cts not present in a bregma but columns added here have NA values
# replace NAs with 0's
gtSpotTopics_N7[is.na(gtSpotTopics_N7)] <- 0
gtSpotTopics_N7[1:3,]

# 3. gtCtGenes
# the `gtCtGenes` is the beta of the average gene expression for each cell cluster,
# in this case using all of cells across the bregma in the animal
dim(gtGexpCellTypes)

# 4. cellCounts
# counts of cells in each simulated spot, but also has spot coordinates for easy plotting
cellCounts_N7 <- do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- simBregmasFN7[[ix]]$cellCounts
  m
}))
dim(cellCounts_N7)

cellCounts_N7[1:10,]

# 5. annotDf
# recall this is meta data data frame includes information for only the cells that were kept in simulated spots
annotDf_N7 <- do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- simBregmasFN7[[ix]]$annotDf
  m
}))
dim(annotDf_N7)

annotDf_N7[1:10,]

# construct the list similar to the output of `buildBregmaCorpus()`
simFN7 <- list(sim = sim_N7,
               gtSpotTopics = gtSpotTopics_N7,
               gtCtGenes = gtGexpCellTypes,
               cellCounts = cellCounts_N7,
               # classColors = classColors,
               annotDf = annotDf_N7)

# num cell types per pixel distribution
table(apply(gtSpotTopics_N7, MARGIN = 1, function(x) length(which(x > 0))))

saveRDS(simFN7, '../UROP/objects/MERFISH_truth_100um2_FULL.rds')
puck <- SpatialRNA(counts = t(as.matrix(sim_N7)), use_fake_coords = T)
puck@coords <- cellCounts_N7[,1:2]
saveRDS(puck, '../UROP/objects/MERFISH_puck_100um2_FULL.rds')

