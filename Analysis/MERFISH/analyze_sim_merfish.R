datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))
source(file.path(datadir, 'R/analysis.R'))

RCTD_list <- readRDS(file.path(datadir, 'Objects/MERFISH_20um2_10_10.rds'))
RCTDpred <- RCTD_list[[length(RCTD_list)]]
results <- readRDS(file.path(datadir, 'Objects/MERFISH_truth_20um2.rds'))

# set up results data frames
pred_results <- RCTDpred@results$results_df
truth_results <- as.data.frame(results$gtSpotTopics[rownames(pred_results),])

doublet_certain <- apply(truth_results, MARGIN = 1, function(x) length(which(x >= 0.2)) == 2)
singlet <- apply(truth_results, MARGIN = 1, function(x) max(x) > 0.8)
first_type <- apply(truth_results, MARGIN = 1, function(x) names(which(x == sort(x, decreasing = TRUE)[1]))[1])
second_type <- apply(truth_results, MARGIN = 1, function(x) names(which(x == sort(x, decreasing = TRUE)[2]))[1])

truth_results$spot_class <- 'reject'
truth_results[doublet_certain,]$spot_class <- 'doublet_certain'
truth_results[singlet,]$spot_class <- 'singlet'
truth_results$first_type <- first_type
truth_results$second_type <- second_type
truth_results <- truth_results[,c('spot_class', 'first_type', 'second_type')]

# reassign doublet/singlet
reassign <- rownames(pred_results[pred_results$spot_class == 'doublet_certain' & 
                                    pred_results$singlet_score - pred_results$min_score < 25,])
pred_results[reassign,]$spot_class <- 'singlet'

# singlets
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

# accuracy
cell_types <- colnames(results$gtSpotTopics)
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


# save results
correct_10_25 <- correct
total_10_25 <- total
doublet_10_25 <- doublet_dist

merfish_accuracy <- data.frame(correct_10_10, correct_10_15, correct_10_20, correct_10_25, 
                                  correct_5_10, correct_5_15, correct_5_20, correct_5_25, 
                                  total_10_10, total_10_15, total_10_20, total_10_25, 
                                  total_5_10, total_5_15, total_5_20, total_5_25)
merfish_doublet <- data.frame(c(doublet_10_10[1:2]), c(doublet_10_15[1:2]), c(doublet_10_20[1:2]), c(doublet_10_25[1:2]),
                                 c(doublet_5_10[1:2]), c(doublet_5_15[1:2]), c(doublet_5_20[1:2]), c(doublet_5_25[1:2]))
colnames(merfish_doublet) <- c('10_10','10_15','10_20','10_25','5_10','5_15','5_20','5_25')

saveRDS(merfish_accuracy,file.path(datadir,'Objects/MERFISH_reassigned_accuracy.rds'))
saveRDS(merfish_doublet,file.path(datadir,'Objects/MERFISH_reassigned_doublet.rds'))


# generate plots

full_accuracy <- readRDS(file.path(datadir, 'Objects/MERFISH_full_accuracy.rds'))
full_doublet <- readRDS(file.path(datadir, 'Objects/MERFISH_full_doublet.rds'))
full_doublet <- full_doublet[c('doublet_certain','singlet'),]
rownames(full_doublet) <- c('singlet','doublet')

# plot doublet distribution
#full_doublet <- apply(X = full_doublet, MARGIN = 2, FUN = function(x) x/sum(x))
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
  ggtitle('Accuracy Across Parameters')


# doublet proportions


