datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))
source(file.path(datadir, 'R/analysis.R'))

RCTD_list <- readRDS(file.path(datadir, 'Objects/cerebellum_10_10.rds'))
RCTDpred <- RCTD_list[[length(RCTD_list)]]

## INITIAL DATA ANALYSIS

# Spatially plot predicted cells
plot_cell_types(RCTDpred)
# Gene expression MSE compare consecutive (testing convergence)
for (i in 1:(length(RCTD_list)-1)) {
  print(gene.mse(RCTD_list[[i]], RCTD_list[[i+1]]))
}
# Pixel assignment accuracy compare consecutive (testing convergence)
for (i in 1:(length(RCTD_list)-1)) {
  print(assignment_accuracy(RCTD_list[[i]], RCTD_list[[i+1]]))
}
# spot type distribution over iterations
for (i in 1:length(RCTD_list)) {
  print(table(RCTD_list[[i]]@results$results_df$spot_class))
}

# compare spot types
RCTD_list <- readRDS(file.path(datadir, 'Objects/cerebellum_10_10.rds'))
RCTDpred <- RCTD_list[[length(RCTD_list)]]
pred_results <- RCTDpred@results$results_df
truth_results <- RCTD_truth@results$results_df
common <- intersect(rownames(pred_results), rownames(truth_results))
mytable <- table(pred_results[common, 'spot_class'], truth_results[common, 'spot_class'])
confmat <- confusionMatrix(pred_results[common, 'spot_class'], truth_results[common, 'spot_class'])
#print(mytable)
#print(confmat$byClass[,'Precision'])


# COMPARED TO RCTD ------------------------

# Accuracy metric for singlets and doublets

RCTD_list <- readRDS(file.path(datadir, 'Objects/cerebellum_10_10.rds'))
RCTD_pred <- RCTD_list[[length(RCTD_list)]]
RCTD_truth <- readRDS(file.path(datadir, 'data/RCTD_cerebellum_slideseq.rds'))
pred_results <- RCTD_pred@results$results_df
truth_results <- RCTD_truth@results$results_df

# merge MLI1 and MLI2 cell types
truth_results$first_type <- gsub('MLI1|MLI2', 'MLI', truth_results$first_type)
truth_results$second_type <- gsub('MLI1|MLI2', 'MLI', truth_results$second_type)

# reassign doublet_like_cutoff
#reassign <- rownames(pred_results[pred_results$spot_class == 'doublet_certain' & 
#                                  pred_results$singlet_score - pred_results$min_score < 25,])
#pred_results[reassign,]$spot_class <- 'singlet'

# assign clusters to ground truth cell types
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
mytable <- singlet_table(truth_results, pred_results)
cluster_assignments <- rownames(mytable)[apply(mytable, 2, function(x) which(x == max(x)))]
names(cluster_assignments) <- colnames(mytable)

# replace cluster with assigned cell type
ground_truth <- truth_results[truth_results$spot_class %in% c('singlet', 'doublet_certain', 'doublet_uncertain'),
                              c('spot_class', 'first_type', 'second_type')]
predicted <- pred_results[rownames(ground_truth),c('spot_class', 'first_type', 'second_type')]
predicted$first_type <- sapply(predicted$first_type, function(x) cluster_assignments[x])
predicted$second_type <- sapply(predicted$second_type, function(x) cluster_assignments[x])

# compute accuracy
cell_types <- RCTD_truth@cell_type_info$info[[2]]
cell_types <- cell_types[! cell_types %in% c('MLI1', 'MLI2')]
cell_types <- append(cell_types, 'MLI')
total <- rep(0, length(cell_types)); correct <- rep(0, length(cell_types))
names(total) <- cell_types; names(correct) <- cell_types

for (barcode in rownames(ground_truth)) {
  truth_row <- ground_truth[barcode,]
  pred_row <- predicted[barcode,]
  if (truth_row$spot_class == 'singlet' || truth_row$spot_class == 'doublet_uncertain') {
    truth_type <- truth_row$first_type
    if (pred_row$spot_class %in% c('singlet', 'doublet_uncertain')) {
      if (truth_type == pred_row$first_type) {
        total[truth_type] <- total[truth_type] + 1
        correct[truth_type] <- correct[truth_type] + 1
      } else {
        total[truth_type] <- total[truth_type] + 1
      }
    } else if (pred_row$spot_class == 'doublet_certain') {
      if (truth_type %in% c(pred_row$first_type, pred_row$second_type)) {
        total[truth_type] <- total[truth_type] + 1
        correct[truth_type] <- correct[truth_type] + 0.5
      } else {
        total[truth_type] <- total[truth_type] + 1
      }
    }
  } else if (truth_row$spot_class == 'doublet_certain') {
    truth_type <- c(truth_row$first_type, truth_row$second_type)
    if (pred_row$spot_class %in% c('singlet', 'doublet_uncertain')) {
      if (pred_row$first_type %in% truth_type) {
        total[pred_row$first_type] <- total[pred_row$first_type] + 1
        correct[pred_row$first_type] <- correct[pred_row$first_type] + 1
      } else {
        total[pred_row$first_type] <- total[pred_row$first_type] + 0.5
        total[pred_row$second_type] <- total[pred_row$second_type] + 0.5
      }
    } else if (pred_row$spot_class == 'doublet_certain') {
      if (pred_row$first_type %in% truth_type) {
        total[pred_row$first_type] <- total[pred_row$first_type] + 1
        correct[pred_row$first_type] <- correct[pred_row$first_type] + 1
      } else {
        total[pred_row$first_type] <- total[pred_row$first_type] + 1
      }
      if (pred_row$second_type %in% truth_type) {
        total[pred_row$second_type] <- total[pred_row$second_type] + 1
        correct[pred_row$second_type] <- correct[pred_row$second_type] + 1
      } else {
        total[pred_row$second_type] <- total[pred_row$second_type] + 1
      }
    }
  }
}

sum(correct) / sum(total) # overall accuracy
correct / total # by class accuracy

# proportion of doublets from true doublets, true singlets

pred_doublets <- predicted[predicted$spot_class == 'doublet_certain',]
doublet_dist <- table(ground_truth[rownames(pred_doublets),]$spot_class)[c('singlet', 'doublet_certain', 'doublet_uncertain')]
barplot(doublet_dist, main = 'True Spot Types of Predicted Doublets', xlab = 'True Spot Type', ylab = 'Count')

reassigned_accuracy <- data.frame(correct_10_10, correct_10_15, correct_10_20, correct_10_25, 
                                  correct_5_10, correct_5_15, correct_5_20, correct_5_25, 
                                  total_10_10, total_10_15, total_10_20, total_10_25, 
                                  total_5_10, total_5_15, total_5_20, total_5_25)
reassigned_doublet <- data.frame(c(doublets_10_10), c(doublets_10_15), c(doublets_10_20), c(doublets_10_25),
                                 c(doublets_5_10), c(doublets_5_15), c(doublets_5_20), c(doublets_5_25))



# compare parameters


full_accuracy <- readRDS(file.path(datadir, 'Objects/cerebellum_reassigned_accuracy.rds'))
full_doublet <- readRDS(file.path(datadir, 'Objects/cerebellum_reassigned_doublet.rds'))

# plot doublet distribution
#full_doublet <- apply(X = full_doublet, MARGIN = 2, FUN = function(x) x/sum(x))
full_doublet <- full_doublet[c('doublet_certain','doublet_uncertain','singlet'),]
full_doublet <- data.matrix(full_doublet)
doublet_data <- melt(full_doublet)
colnames(doublet_data) <- c('spot_type', 'params', 'count')
ggplot(doublet_data, aes(fill=spot_type, y=count, x=params)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab('Parameters (QL_score_cutoff, doublet_like_cutoff)') + ylab('Count') +
  ggtitle('True Spot Types of Predicted Doublets') + labs(fill = 'True Spot Type')

# plot overall accuracy 
totals <- colSums(full_accuracy)
for (i in 1:8) {
  totals[i] <- totals[i] / totals[i+8]
}
totals <- totals[1:8]
total_data <- data.frame(totals)
total_data$QL_score_cutoff <- c('10','10','10','10','5','5','5','5')
total_data$doublet_like_cutoff <- c('10','15','20','25','10','15','20','25')

ggplot(data=total_data, aes(x=doublet_like_cutoff, y=totals, group=QL_score_cutoff)) +
  geom_line(aes(color = QL_score_cutoff))+
  geom_point() +
  ylab('Accuracy') +
  ggtitle('Accuracy Across Parameters') +
  ylim(c(0.880, NA))

# plot per cell type accuracy
per_cell <- full_accuracy
for (i in 1:8) {
  per_cell[,i] <- per_cell[,i] / per_cell[,i+8]
}
per_cell <- per_cell[,1:8]
cell_data <- melt(data.matrix(per_cell))
cell_data$QL_score_cutoff <- c(rep('10', 72), rep('5', 72))
cell_data$doublet_like_cutoff <- rep(2:5 * 5, each = 18, times = 2)
cell_data$Var2 <- gsub('correct_', '', cell_data$Var2)

ggplot(data=cell_data[cell_data$QL_score_cutoff == 5,], aes(x=doublet_like_cutoff, y=value, group=Var1, color=Var1)) +
  geom_jitter(width = 0.3, size = 1) +
  ylab('Accuracy') +
  ggtitle('Accuracy By Cell Type (QL_score_cutoff = 5)') +
  labs(fill = 'Cell Type') +
  ylim(c(0.5,1))

# plot bergmann/purkinje marker genes (obtained from marker_gene_analysis.R)
berg_data <- data.frame(c(1.393259,1.466302,1.63852,2.125758,1.466026,1.520642,1.652017,2.031301),
                        c(rep('10', 4), rep('5', 4)),
                        rep(5:2 * 5, times = 2))
colnames(berg_data) <- c('value', 'QL_score_cutoff', 'doublet_like_cutoff')
purk_data <- data.frame(c(1.740934,1.769885,1.922802,2.181597,1.706764,1.730996,1.831023,1.949197),
                        c(rep('10', 4), rep('5', 4)),
                        rep(5:2 * 5, times = 2))
colnames(purk_data) <- c('value', 'QL_score_cutoff', 'doublet_like_cutoff')

ggplot(data=purk_data, aes(x=doublet_like_cutoff, y=value, group=QL_score_cutoff)) +
  geom_line(aes(color = QL_score_cutoff))+
  geom_point() +
  ylab('Marker Gene Expression Log Difference') +
  ggtitle('Purkinje Marker Gene Expression')

# FIND MARKER GENES
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

RCTD_pred <- readRDS(file.path(datadir, 'Objects/cerebellum_10_10.rds'))[[17]]
RCTD_truth <- readRDS(file.path(datadir, 'data/RCTD_cerebellum_slideseq.rds'))
pred_results <- RCTD_pred@results$results_df
truth_results <- RCTD_truth@results$results_df
# merge MLI1 and MLI2 cell types
truth_results$first_type <- gsub('MLI1|MLI2', 'MLI', truth_results$first_type)
truth_results$second_type <- gsub('MLI1|MLI2', 'MLI', truth_results$second_type)
# assign clusters to ground truth cell types
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
mytable <- singlet_table(truth_results, pred_results)
cluster_assignments <- rownames(mytable)[apply(mytable, 2, function(x) which(x == max(x)))]
names(cluster_assignments) <- colnames(mytable)

myRCTD <- readRDS(file.path(datadir, 'Objects/cerebellum_10_10.rds'))[[17]]
cur_cell_types <- c("2","0","4","3","5") # set to cell types
puck <- myRCTD@spatialRNA
cell_type_info_restr = myRCTD@cell_type_info$info
cell_type_info_restr[[1]] = cell_type_info_restr[[1]][,cur_cell_types]
cell_type_info_restr[[2]] = cur_cell_types; cell_type_info_restr[[3]] = length(cur_cell_types)
de_genes <- get_de_genes(cell_type_info_restr, puck, fc_thresh = 3, expr_thresh = .0001, MIN_OBS = 3)
marker_data_de = get_marker_data(cell_type_info_restr[[2]], cell_type_info_restr[[1]], de_genes)

berg_pred <- rownames(marker_data_de)[marker_data_de$cell_type == "2"]
purk_pred <- rownames(marker_data_de)[marker_data_de$cell_type == "5"]

marker_data_truth <- readRDS(file.path(datadir, 'Data/marker_data_de_standard.RDS'))
berg_truth <- rownames(marker_data_truth)[marker_data_truth$cell_type == "Bergmann"]
purk_truth <- rownames(marker_data_truth)[marker_data_truth$cell_type == "Purkinje"]
