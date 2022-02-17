datadir <- '~/UROP'


# GENERATE SIMULATED PUCK WITH 1000 UMI, UNIFORM CELL TYPE MIXING ---------------------------------

sub_sample_cell <- function(gene_list, raw.data, cell_index, UMI_sample) {
  sub_sample = sample(rep(rownames(raw.data),raw.data[,cell_index]),UMI_sample,replace=FALSE)
  sub_sample = table(sub_sample)[as.character(gene_list)]
  names(sub_sample) = gene_list
  sub_sample[is.na(sub_sample)] <- 0
  return(sub_sample)
}

generate_sim_puck <- function(common_cell_types, gene_list, ref, trials = 30, UMI_tot = 1000) {
  n_cell_types = length(common_cell_types)
  #trials = 30
  n_conditions = 13
  boundary = ceiling(n_conditions / 2) # DE occurs from conditions 7-13
  N_samples = (n_cell_types * trials * n_conditions * (n_cell_types - 1))/2 # total number of samples
  first_UMI = numeric(N_samples); first_type = character(N_samples); second_type = character(N_samples)
  # first_UMI, first_type, second_type are all vectors of length N_samples. first_type, second_type store cell types,
  # first_UMI stores UMI counts
  UMI_tot = UMI_tot; UMI_step = UMI_tot / (n_conditions-1)
  UMI_tot = round(UMI_step * (n_conditions-1)); UMI1_vec = round(0:(n_conditions-1)*UMI_step)
  # above are the UMI amounts
  beads = Matrix(0, nrow = length(gene_list), ncol = N_samples)
  rownames(beads) = gene_list; colnames(beads) = 1:N_samples
  # beads that are filled as the algorithm runs, only thing that has to change is N_samples
  firstbeads = Matrix(0, nrow = length(gene_list), ncol = N_samples)
  rownames(firstbeads) = gene_list; colnames(firstbeads) = 1:N_samples
  secbeads = Matrix(0, nrow = length(gene_list), ncol = N_samples)
  rownames(secbeads) = gene_list; colnames(secbeads) = 1:N_samples
  # firstbeads and secbeads add together to create beads
  first_index_list <- numeric(N_samples); second_index_list <- numeric(N_samples);
  names(first_index_list) <- 1:N_samples; names(second_index_list) <- 1:N_samples;
  index = 1
  nUMI = ref@nUMI
  results <- matrix(ncol = 4)
  for(i in 1:(n_cell_types-1)) {
    print(paste("Progress",i))
    for(j in (i+1):n_cell_types) {
      print(paste("ProgressSecond",j))
      type1 = common_cell_types[i]; type2 = common_cell_types[j]
      for (condition in 1:n_conditions) {
        UMI1 = UMI1_vec[condition]; UMI2 = UMI_tot - UMI1
        for(t in 1:trials) {
          first_UMI[index] = UMI1; first_type[index] = type1; second_type[index] = type2
          # sets UMI values and cell types
          firstInd = sample(intersect(which(ref@cell_types == type1) , which(nUMI > 1000)),1)
          secondInd = sample(intersect(which(ref@cell_types == type2) , which(nUMI > 1000)),1)
          # samples a reference barcode index for each cell type
          first_index_list[index] <- firstInd; second_index_list[index] <- secondInd
          # stores the barcode indices in the lists
          firstbeads[,index] = as.vector(sub_sample_cell(gene_list, ref@counts, firstInd, UMI1))
          secbeads[,index] = as.vector(sub_sample_cell(gene_list, ref@counts, secondInd, UMI2))
          beads[,index] = firstbeads[,index] + secbeads[,index]
          # samples gene counts using the UMI amounts and the barcode index from above
          results <- rbind(results, c(type1, type2, UMI1, UMI2))
          index = index + 1
        }
      }
    }
  }
  #UMI_vect <- rep(UMI_tot,N_samples)
  #names(UMI_vect) <- colnames(beads)
  #puck <- SpatialRNA(NULL, beads, nUMI = UMI_vect, use_fake_coords = T)
  results <- results[-1,]
  results <- as.data.frame(results)
  return(list(beads = beads, results = results))
}

reference <- readRDS(file.path(datadir, '/data/reference_RCTD_vec.rds'))
cell_types <- c('Astrocytes', 'Bergmann', 'Fibroblast', 'Golgi', 'Granule', 'MLI1', 'MLI2', 'Oligodendrocytes', 'Polydendrocytes', 'Purkinje')
puck <- generate_sim_puck(cell_types, rownames(reference@counts), reference, trials=30)
saveRDS(puck, file.path(datadir, '/data/sim_puck.rds'))

# generate ground truth types (now included in generate_sim_puck)
puck <- readRDS(file.path(datadir, '/data/sim_puck.rds'))
n_cell_types = length(cell_types)
trials = 30
n_conditions = 13
boundary = ceiling(n_conditions / 2)
N_samples = (n_cell_types * trials * n_conditions * (n_cell_types - 1))/2
first_UMI = numeric(N_samples); first_type = character(N_samples); second_type = character(N_samples)
UMI_tot = 1000; UMI_step = UMI_tot / (n_conditions-1)
UMI_tot = round(UMI_step * (n_conditions-1)); UMI1_vec = round(0:(n_conditions-1)*UMI_step)
barcodes <- colnames(puck@counts)
results <- matrix(ncol = 4)
colnames(results) <- c('type1', 'type2', 'UMI1', 'UMI2')
for(i in 1:(n_cell_types-1)) {
  for(j in (i+1):n_cell_types) {
    type1 = cell_types[i]; type2 = cell_types[j]
    for (condition in 1:n_conditions) {
      UMI1 = UMI1_vec[condition]; UMI2 = UMI_tot - UMI1
      for(t in 1:trials) {
        results <- rbind(results, c(type1, type2, UMI1, UMI2))
      }
    }
  }
}
results <- results[-1,]
results <- as.data.frame(results)
saveRDS(results, file.path(datadir, '/data/sim_puck_truth.rds'))


# GENERATE SIMULATED PUCK WITH 1000 UMI, PAIRS OF CELL TYPES MIXING ---------------------------------
# Cell types mixed: Purkinje/Bergmann, MLI1/MLI2/Granule, Astrocytes/Oligodendrocytes

reference <- readRDS(file.path(datadir, '/data/reference_RCTD_vec.rds'))
cell_groups <- list(c('Bergmann', 'Purkinje'), c('Granule', 'MLI1', 'MLI2'), c('Astrocytes', 'Oligodendrocytes'))
for (i in 1:length(cell_groups)) {
  cell_types <- cell_groups[[i]]
  beads_list <- generate_sim_puck(cell_types, rownames(reference@counts), reference, trials=30)
  if (i == 1) {
    beads <- beads_list$beads
    results <- beads_list$results
  } else {
    beads <- cbind(beads, beads_list$beads)
    results <- rbind(results, beads_list$results)
  }
}
colnames(beads) <- 1:dim(beads)[2]
rownames(results) <- 1:dim(beads)[2]
puck <- SpatialRNA(NULL, beads, use_fake_coords = T)
saveRDS(puck, file.path(datadir, '/objects/sim_puck_paired.rds'))
saveRDS(results, file.path(datadir, '/objects/sim_puck_paired_truth.rds'))


# GENERATE SIMULATED PUCK WITH VARIABLE nUMI, UNIFORM MIXING ---------------------------------

reference <- readRDS(file.path(datadir, '/data/reference_RCTD_vec.rds'))
cell_types <- c('Astrocytes', 'Bergmann', 'Fibroblast', 'Golgi', 'Granule', 'MLI1', 'MLI2', 'Oligodendrocytes', 'Polydendrocytes', 'Purkinje')
for (i in 1:10) {
  UMI_tot = 100 * i
  beads_list <- generate_sim_puck(cell_types, rownames(reference@counts), reference, trials=10, UMI_tot=UMI_tot)
  if (i == 1) {
    beads <- beads_list$beads
    results <- beads_list$results
  } else {
    beads <- cbind(beads, beads_list$beads)
    results <- rbind(results, beads_list$results)
  }
}
colnames(beads) <- 1:dim(beads)[2]
rownames(results) <- 1:dim(beads)[2]
puck <- SpatialRNA(NULL, beads, use_fake_coords = T)
saveRDS(puck, file.path(datadir, '/objects/sim_puck_variableUMI.rds'))
saveRDS(results, file.path(datadir, '/objects/sim_puck_variableUMI_truth.rds'))

