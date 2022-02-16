# note to self: this classification metric works poorly with small numbers of clusters, since different cell types will necessarily be grouped together

# cleaning tables
mytable <- cell.table(RCTDtruth, RCTDpred)
mytable <- rbind(mytable, mytable['MLI1',] + mytable['MLI2',])
rownames(mytable)[length(rownames(mytable))] <- 'MLI'
mytable <- mytable[-c(which(rownames(mytable) == 'MLI1'), which(rownames(mytable) == 'MLI2')),]

# new accuracy metric
assignment.accuracy <- function(mytable, N) {
	cluster_assignments <- apply(mytable, 2, function(x) which(x == max(x)))
	truth_totals <- c(dim(mytable)[1])
	for (i in 1:dim(mytable)[1]) {
		truth_totals[i] <- sum(mytable[i,] * (cluster_assignments == i))
	}
	cell_accuracy <- truth_totals / apply(mytable, 1, sum)
	cell_accuracy <- sort(cell_accuracy, decreasing=TRUE)
	return(mean(cell_accuracy[1:N]))
}


gene_expression <- gene.expression(RCTDtruth, RCTDlist[[6]])
mytable <- cell.table(RCTDtruth, RCTDlist[[6]])
cluster_assignments <- apply(mytable, 2, function(x) which(x == max(x)))
pred_ranks = c()
ref_ranks = c()
for (i in 1:length(cluster_assignments)) {
	pred_cell = i
	truth_cel = cluster_assignments[i]
	rank <- gene.expression.rank(gene_expression, truth_cell, pred_cell)
	pred_ranks <- append(pred_ranks, rank$prediction)
	ref_ranks <- append(ref_ranks, rank$reference)
}

plot(ref_expression, pred_expression, xlab="Reference Gene Expression Rank", ylab="Predicted Gene Expression Rank")

gene.expression.rank <- function(gene_expression, truth_cell, pred_cell)
	ref_expression <- sort(gene_expression$reference[,truth_cell], decreasing = TRUE)
	pred_expression <- sort(gene_expression$prediction[,pred_cell], decreasing = TRUE)
	for (i in 1:length(ref_expression)) {
	    ref_expression[i] = i
	}
	for (i in 1:length(pred_expression)) {
	    pred_expression[i] = i
	}
	ref_expression <- ref_expression[order(factor(names(ref_expression)))]
	pred_expression <- pred_expression[order(factor(names(pred_expression)))]
	return(list(reference = ref_expression, prediction = pred_expression))

plot(ref_expression, pred_expression)
plot(ref_expression, pred_expression, xlim = c(0,100), ylim = c(0,100))



# compare gene expression of stdeconvolve with ground truth with iterative alg
results <- readRDS('../UROP/objects/stdeconvolve_cerebellum_slideseq_23res.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
ref_info <- RCTDtruth@cell_type_info$renorm[[1]]
pred_info <- t(results$beta)
gene_list <- intersect(rownames(ref_info), rownames(pred_info))
ref_info <- data.matrix(ref_info[gene_list,])
pred_info <- data.matrix(pred_info[gene_list,])

aggregate.gene.rank


# compare marker gene expression
marker_data_de <- readRDS('../UROP/data/marker_data_de_standard.RDS')
berg_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Bergmann"]
purk_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Purkinje"]
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
RCTDlist <- readRDS('../UROP/objects/RCTD_list_truth.rds')
RCTDpred <- RCTDlist[[1]]
cell.table(RCTDtruth, RCTDpred)

berg_log_df <- log(RCTDpred@cell_type_info$renorm[[1]][berg_genes,c(3,6)])
purk_log_df <- log(RCTDpred@cell_type_info$renorm[[1]][purk_genes,c(3,6)])
colnames(berg_log_df) <- c('Bergmann', 'Purkinje')
colnames(purk_log_df) <- c('Bergmann', 'Purkinje')
x <- berg_log_df$Bergmann-berg_log_df$Purkinje
x <- x[!is.na(x) & !is.infinite(x)]
mean(x)
x <- purk_log_df$Purkinje-purk_log_df$Bergmann
x <- x[!is.na(x) & !is.infinite(x)]
mean(x)
berg_log_df$type <- 'Bergmann'
purk_log_df$type <- 'Purkinje'
full_df <- rbind(berg_log_df, purk_log_df)

my_pal = pals::coolwarm(20)
p1 <- ggplot2::ggplot(full_df) + 
    ggplot2::geom_point(ggplot2::aes(x=Bergmann,y=Purkinje, color = type),alpha = 0.2,size=1) + 
    ggplot2::coord_fixed() + ggplot2::theme_classic()  +
    labs(color="Marker Cell Type") + guides(color = guide_legend(override.aes = list(size = 3,alpha=1))) + 
    scale_color_manual(values=c(my_pal[1], my_pal[20]))+ theme(legend.position="top") + xlab('Bergmann Log Expression') + ylab('Purkinje Log Expression') +
    geom_abline(intercept = 0, slope = 1, color="black", size=0.5, linetype="dashed")
p1

#same but for ground truth
marker_data_de <- readRDS('../UROP/data/marker_data_de_standard.RDS')
berg_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Bergmann"]
purk_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Purkinje"]
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')

berg_log_df <- log(RCTDtruth@cell_type_info$renorm[[1]][berg_genes,c('Bergmann','Purkinje')])
purk_log_df <- log(RCTDtruth@cell_type_info$renorm[[1]][purk_genes,c('Bergmann','Purkinje')])
x <- berg_log_df$Bergmann-berg_log_df$Purkinje
x <- x[!is.na(x) & !is.infinite(x)]
mean(x)
x <- purk_log_df$Purkinje-purk_log_df$Bergmann
x <- x[!is.na(x) & !is.infinite(x)]
mean(x)
berg_log_df$type <- 'Bergmann'
purk_log_df$type <- 'Purkinje'
full_df <- rbind(berg_log_df, purk_log_df)

my_pal = pals::coolwarm(20)
p1 <- ggplot2::ggplot(full_df) + 
    ggplot2::geom_point(ggplot2::aes(x=Bergmann,y=Purkinje, color = type),alpha = 0.2,size=1) + 
    ggplot2::coord_fixed() + ggplot2::theme_classic()  +
    labs(color="Marker Cell Type") + guides(color = guide_legend(override.aes = list(size = 3,alpha=1))) + 
    scale_color_manual(values=c(my_pal[1], my_pal[20]))+ theme(legend.position="top") + xlab('Bergmann Log Expression') + ylab('Purkinje Log Expression') +
    geom_abline(intercept = 0, slope = 1, color="black", size=0.5, linetype="dashed")
p1

#same but for stdeconvolve results
marker_data_de <- readRDS('../UROP/data/marker_data_de_standard.RDS')
berg_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Bergmann"]
purk_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Purkinje"]
results <- readRDS('../UROP/objects/stdeconvolve_cerebellum_slideseq_13res.rds')
RCTDtruth <- readRDS('../UROP/objects/RCTD_cerebellum_slideseq.rds')
berg_log_df <- log(t(results$beta)[intersect(rownames(results$beta),berg_genes),c(4,1)])
purk_log_df <- log(t(results$beta)[intersect(rownames(results$beta),purk_genes),c(4,1)])
colnames(berg_log_df) <- c('Bergmann', 'Purkinje')
colnames(purk_log_df) <- c('Bergmann', 'Purkinje')
berg_log_df$type <- 'Bergmann'
purk_log_df$type <- 'Purkinje'
full_df <- rbind(berg_log_df, purk_log_df)

my_pal = pals::coolwarm(20)
p1 <- ggplot2::ggplot(full_df) + 
    ggplot2::geom_point(ggplot2::aes(x=Bergmann,y=Purkinje, color = type),alpha = 0.2,size=1) + 
    ggplot2::coord_fixed() + ggplot2::theme_classic()  +
    labs(color="Marker Cell Type") + guides(color = guide_legend(override.aes = list(size = 3,alpha=1))) + 
    scale_color_manual(values=c(my_pal[1], my_pal[20]))+ theme(legend.position="top") + xlab('Bergmann Log Expression') + ylab('Purkinje Log Expression') +
    geom_abline(intercept = 0, slope = 1, color="black", size=0.5, linetype="dashed")
p1
