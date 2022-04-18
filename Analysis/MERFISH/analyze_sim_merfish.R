datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))
source(file.path(datadir, 'R/analysis.R'))

RCTD_list <- readRDS(file.path(datadir, 'Objects/MERFISH_20um2_10_5.rds'))
RCTDpred <- RCTD_list[[length(RCTD_list)]]

results <- readRDS(file.path(datadir, 'Objects/MERFISH_truth_20um2.rds'))

RCTDpred <- readRDS(file.path(datadir, 'Objects/MERFISH_20um2_10_5.rds'))

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
#reassign <- rownames(pred_results[pred_results$spot_class == 'doublet_certain' & 
#                                    pred_results$singlet_score - pred_results$min_score < 25,])
#pred_results[reassign,]$spot_class <- 'singlet'

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
correct_sil_80 <- correct
total_sil_80 <- total
doublet_sil_80 <- doublet_dist

merfish_accuracy <- data.frame(correct_sil_00, correct_sil_20, correct_sil_40, correct_sil_60, correct_sil_80,
                               total_sil_00, total_sil_20, total_sil_40, total_sil_60, total_sil_80)
merfish_doublet <- data.frame(c(doublet_sil_00[1:2]), c(doublet_sil_20[1:2]), c(doublet_sil_40[1:2]), c(doublet_sil_60[1:2]),
                                 c(doublet_sil_80[1:2]))
colnames(merfish_doublet) <- c('0.00','0.20','0.40','0.60','0.80')

saveRDS(merfish_accuracy,file.path(datadir,'Objects/MERFISH_silhouette_accuracy.rds'))
saveRDS(merfish_doublet,file.path(datadir,'Objects/MERFISH_silhouette_doublet.rds'))


# generate plots

full_accuracy <- readRDS(file.path(datadir, 'Objects/MERFISH_silhouette_accuracy.rds'))
full_doublet <- readRDS(file.path(datadir, 'Objects/MERFISH_silhouette_doublet.rds'))
full_doublet <- full_doublet[c('doublet_certain','singlet'),]
rownames(full_doublet) <- c('singlet','doublet')

# plot doublet distribution
#full_doublet <- apply(X = full_doublet, MARGIN = 2, FUN = function(x) x/sum(x))
full_doublet <- data.matrix(full_doublet)
doublet_data <- melt(full_doublet)
colnames(doublet_data) <- c('spot_type', 'params', 'count')
doublet_data$params <- as.character(doublet_data$params)
ggplot(doublet_data, aes(fill=spot_type, y=count, x=params)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab('Silhouette Cutoff') + ylab('Count') +
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





# FULL MODE COMPARISON
RCTD_list <- readRDS(file.path(datadir, 'Objects/MERFISH_100um2_fullmode.rds'))
RCTDpred <- RCTD_list[[100]]
results <- readRDS(file.path(datadir, 'Objects/MERFISH_truth_100um2.rds'))

pred_results <- RCTDpred@results$weights
for (i in 1:dim(pred_results)[1]) {
  pred_results[i,] = pred_results[i,] / rowSums(pred_results)[i]
}
truth_results <- results$gtSpotTopics[rownames(pred_results),]

#  analyze majority cell type in each pixel
pred_maj <- colnames(pred_results)[apply(pred_results, 1, which.max)]
truth_maj <- colnames(truth_results)[apply(truth_results, 1, which.max)]
singlet_table <- table(pred_maj, truth_maj) # does not work very well...

# analyze gene expression
pred_gene <- RCTDpred@cell_type_info$renorm[[1]]
truth_gene <- t(results$gtCtGenes)
for (i in 1:dim(pred_gene)[2]) {
  pred_gene[,i] = pred_gene[,i] / colSums(pred_gene)[i]
}

pred_gene <- log(pred_gene)
truth_gene <- log(truth_gene)
# correlation
expression <- as.matrix(cbind(pred_gene, truth_gene))
expression <- expression[rowSums(is.infinite(expression)) == 0, ]
correlation <- cor(expression)[11:18,1:10]
data <- melt(correlation ^ 2)
ggplot(data, aes(factor(Var1), factor(Var2), fill= value)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colors = pals::brewer.blues(20)[2:20], limits=c(0,1), name='r^2') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Reference Cell Type')+ ylab('Predicted Cell Type')
# MSE and KL divergence
ref_expression <- data.matrix(truth_gene)
pred_expression <- data.matrix(pred_gene)
gene_list <- intersect(rownames(na.omit(ref_expression[rowSums(is.infinite(ref_expression)) == 0,])), rownames(na.omit(pred_expression[rowSums(is.infinite(pred_expression)) == 0,])))
ref_expression <- ref_expression[gene_list,]
pred_expression <- pred_expression[gene_list,]
ref_cells <- colnames(ref_expression); pred_cells <- colnames(pred_expression)
mse_mat <- matrix(,nrow = length(ref_cells), ncol = length(pred_cells))
rownames(mse_mat) <- ref_cells; colnames(mse_mat) <- pred_cells
kl_mat <- matrix(,nrow = length(ref_cells), ncol = length(pred_cells))
rownames(kl_mat) <- ref_cells; colnames(kl_mat) <- pred_cells
for (i in 1:length(ref_cells)) {
  for (j in 1:length(pred_cells)) {
    mse_mat[i,j] <- sum((ref_expression[,i] - pred_expression[,j]) ** 2) / length(gene_list)
    kl_mat[i,j] <- sum((ref_expression[,i] - pred_expression[,j]) * exp(ref_expression[,i])) / length(gene_list)
  }
}

majority_assign <- colnames(singlet_table)[apply(singlet_table, 1, function(x) which(x == max(x)))]
corr_assign <- rownames(correlation)[apply(correlation, 2, function(x) which(x == max(x)))]
mse_assign <- rownames(mse_mat)[apply(mse_mat, 2, function(x) which(x == min(x)))]
kl_assign <- rownames(kl_mat)[apply(kl_mat, 2, function(x) which(x == min(x)))]
majority_assign
corr_assign
mse_assign
kl_assign

# plot puck
plot_puck_continuous(
  RCTDpred@originalSpatialRNA,
  colnames(RCTDpred@originalSpatialRNA@counts),
  RCTDpred@results$weights[,6],
  size=5,
  my_pal = pals::brewer.blues(20)[2:20]
)

# corelation matrix and heatmap
cor(cbind(as.matrix(pred_results), as.matrix(truth_results)))[1:0,11:18]
data <- melt(abs(correlation))
ggplot(data, aes(factor(Var1), factor(Var2), fill= value)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colors = pals::brewer.blues(20)[2:20], limits=c(0,1), name='r^2') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Reference Cell Type')+ ylab('Predicted Cell Type')

# scatter plots
data <- as.data.frame(cbind(pred_results[,1] + pred_results[,2] + pred_results[,3] + pred_results[,4], truth_results[,'Inhibitory']))
colnames(data) <- c('pred','truth')

my_pal = pals::coolwarm(20)
p2 <- ggplot(data, aes(x=truth, y=pred, colour="blue")) +
  geom_point(color=my_pal[1]) +
  theme_classic() +
  xlab('True Inhibitory Proportion') +
  ylab('Predicted Inhibitory Proportion') + 
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-.03,1.03)) + 
  scale_x_continuous(breaks = c(0,0.5,1), limits = c(-.03,1.03)) +
  theme(legend.position = "none") +
  geom_abline(slope=1, color = my_pal[20])
p2

sqrt(sum((pred_results[,1] + pred_results[,2] + pred_results[,3] + pred_results[,4] - truth_results[,'Inhibitory']) ** 2) / length(pred_results[,8]))
