datadir <- '~/UROP'
source(file.path(datadir, 'R/algorithm.R'))
source(file.path(datadir, 'R/analysis.R'))


# load data
puck <- readRDS(file.path(datadir, 'objects/sim_puck_paired.rds'))
results <- readRDS(file.path(datadir, 'objects/sim_puck_paired_truth.rds'))
colnames(results) <- c('type1','type2','UMI1','UMI2')
RCTD_list <- readRDS(file.path(datadir, 'objects/sim_puck_paired_results_8res_SCT_SILHOUETTE.rds'))
RCTDpred <- RCTD_list[[length(RCTD_list)]]
pred_results <- RCTDpred@results$results_df

# analyze singlets
truth_types = c()
pred_types = c()
for (i in 1:dim(results)[1]) {
  if (pred_results[i,'spot_class'] == 'singlet' & strtoi(results[i,'UMI1']) > strtoi(results[i,'UMI2'])) {
    truth_types <- append(truth_types, results[i,'type1'])
    pred_types <- append(pred_types, RCTDpred@results$results_df[i, 'first_type'])
  } else if (pred_results[i,'spot_class'] == 'singlet' & strtoi(results[i,'UMI2']) > strtoi(results[i,'UMI1'])) {
    truth_types <- append(truth_types, results[i,'type2'])
    pred_types <- append(pred_types, RCTDpred@results$results_df[i, 'first_type'])
  }
}
#truth_types <- factor(truth_types, levels = cell_types)
#pred_types <- factor(pred_types, levels = levels(RCTDpred@results$results_df$first_type))
mytable <- table(truth_types, pred_types)
mytable
assignment.accuracy(mytable)
assignment.accuracy.plot(mytable)

# PLOT DOUBLETS FOR BERG/PURK

berg_purk_truth <- results[results$type1 == 'Astrocytes' & results$type2 == 'Oligodendrocytes',]
berg_purk_pred <- pred_results[rownames(berg_purk_truth),]

#berg_purk_pred_singlets <- berg_purk_pred[berg_purk_pred$spot_class == 'singlet' & (berg_purk_pred$first_type == '4' | berg_purk_pred$first_type == '8'),]
berg_purk_pred_doublets <- berg_purk_pred[(berg_purk_pred$spot_class == 'doublet_certain') &
                 ((berg_purk_pred$first_type == '4' & berg_purk_pred$second_type == '5') | (berg_purk_pred$first_type == '5' & berg_purk_pred$second_type == '4')),]
berg_first_IDs <- rownames(berg_purk_pred_doublets[berg_purk_pred_doublets$first_type == '5',])
berg_second_IDs <- rownames(berg_purk_pred_doublets[berg_purk_pred_doublets$second_type == '5',])

berg_first_IDs <- rownames(berg_purk_pred[berg_purk_pred$first_type == '5',])
berg_second_IDs <- rownames(berg_purk_pred[berg_purk_pred$second_type == '5',])

berg_first_truth <- strtoi(berg_purk_truth[berg_first_IDs, 'UMI1']) / 1000
berg_second_truth <- strtoi(berg_purk_truth[berg_second_IDs, 'UMI1']) / 1000
berg_first_pred <- RCTDpred@results$weights_doublet[berg_first_IDs, 'first_type']
berg_second_pred <- RCTDpred@results$weights_doublet[berg_second_IDs, 'second_type']

berg_true_prop <- c(berg_first_truth, berg_second_truth)
berg_pred_prop <- c(berg_first_pred, berg_second_pred)


# MEAN SQUARED ERROR OF DOUBLET PROPORTIONS, DOUBLET CLASSIFICATION RATE

sqrt(sum((berg_true_prop - berg_pred_prop)^2) / length(berg_true_prop))

plot(berg_true_prop, berg_pred_prop, xlim = c(0,1), ylim = c(0,1), xlab = 'True Bergmann proportion', ylab = 'Predicted Bergmann proportion')

prop_df <- data.frame(factor(berg_true_prop), berg_pred_prop)
true_proportions <- levels(prop_df$factor.berg_true_prop.)
means <- c()
stdevs <- c()
for (prop in true_proportions) {
  means <- append(means, mean(prop_df[prop_df$factor.berg_true_prop. == prop,]$berg_pred_prop))
  stdevs <- append(stdevs, sd(prop_df[prop_df$factor.berg_true_prop. == prop,]$berg_pred_prop))
}
library(ggplot2)
results_df <- data.frame(as.numeric(true_proportions), means, stdevs)
colnames(results_df) <- c('true', 'means', 'stdev')
my_pal = pals::coolwarm(20)

p1 <-ggplot2::ggplot(results_df, ggplot2::aes(x=true, y=means, colour = "true")) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::geom_line(ggplot2::aes(y=true,colour = "means")) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=means-1.96*stdevs, ymax=means+1.96*stdevs), width=.05,
                         position=ggplot2::position_dodge(0.05)) + theme_classic() + xlab('True Bergmann Proportion')+ ylab('Predicted Bergmann Proportion')+ scale_color_manual(values=c(my_pal[20], my_pal[1]),labels = c("",""), name = "") + scale_y_continuous(breaks = c(0,0.5,1), limits = c(-.03,1.03))+ scale_x_continuous(breaks = c(0,0.5,1), limits = c(-.03,1.03))+ theme(legend.position = "none")
p1


# CLASSIFICATION PERCENTAGES

dim(berg_purk_pred[berg_purk_pred$spot_class == 'singlet',])
dim(berg_purk_pred[berg_purk_pred$spot_class == 'singlet' & 
               (berg_purk_pred$first_type == '5' | berg_purk_pred$first_type == '8'),])
dim(berg_purk_pred[berg_purk_pred$spot_class == 'doublet_certain',])
dim(berg_purk_pred[berg_purk_pred$spot_class == 'doublet_certain' &
                   ((berg_purk_pred$first_type == '5' & berg_purk_pred$second_type == '8') | (berg_purk_pred$first_type == '8' & berg_purk_pred$second_type == '5')),])
dim(berg_purk_pred[berg_purk_pred$spot_class == 'doublet_uncertain',])
dim(berg_purk_pred[berg_purk_pred$spot_class == 'doublet_uncertain' &
                     ((berg_purk_pred$first_type == '5' & berg_purk_pred$second_type == '8') | (berg_purk_pred$first_type == '8' & berg_purk_pred$second_type == '5')),])


# DOUBLET CLASSIFICATION RATE

#RCTD_list <- readRDS('../UROP/results/simpuck_results/sim_puck_results_de_fit.rds')
#RCTDpred <- RCTD_list[[7]]
#pred_results <- RCTDpred@results$results_df

doublets_pred <- pred_results[pred_results$spot_class == 'doublet_certain',] #results filtered for confident doublets
doublets_true <- results[rownames(doublets_pred),] #real proportions of predicted doublets
results$prop <- factor(strtoi(pmin(results$UMI1, results$UMI2)) / 1000) #real proportinos of all doublets
doublets_true$prop <- factor(strtoi(pmin(doublets_true$UMI1, doublets_true$UMI2)) / 1000, levels=levels(results$prop))
doublet_classification_rate <- table(doublets_true$prop) / table(results$prop)

doublet_df <- as.data.frame(doublet_classification_rate)
doublet_df$std <- sqrt(doublet_classification_rate * (1 - doublet_classification_rate) / table(results$prop))
colnames(doublet_df) <- c('prop', 'value', 'std')
doublet_df$prop <- as.numeric(levels(doublet_df$prop))

my_pal = pals::coolwarm(20)
p2 <- ggplot(doublet_df, aes(prop,value)) + geom_line()  + 
  geom_errorbar(aes(ymin=value-1.96*std, ymax=value+1.96*std), width=.02,position=position_dodge(.001)) + theme_classic() + 
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.5) + xlab("UMI Proportion of Minority Cell Type") + ylab("Doublet Classification Rate") + scale_color_manual(values=c(my_pal[1], my_pal[20]),labels = c("Doublet"), name = "Class") + scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1))+ scale_x_continuous(breaks = c(0,0.25,0.5), limits = c(-.03,0.53))
p2

# compare doublet rate across all generation methods
#df1$group <- 'CSIDE with singlets and doublets'
#df2$group <- 'CSIDE with singlets'
#df3$group <- 'Mean of singlets'
#df <- rbind(df1, df2, df3)
#p <- ggplot(df, aes(x=prop, y=value, group=group, col=group)) + geom_line()  + 
#  geom_errorbar(aes(ymin=value-1.96*std, ymax=value+1.96*std), width=.02,position=position_dodge(.001)) + theme_classic() + labs(col = "Gene Expression Profile Generation Method") + 
#  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.5) + xlab("UMI Proportion of Minority Cell Type") + ylab("Doublet Classification Rate") + scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1))+ scale_x_continuous(breaks = c(0,0.25,0.5), limits = c(-.03,0.53))
#p


## SINGLET SCORE ANAYLSES

## singlet score vs doublet proportion
# singlet score vs. purkinje proportion for singlets
purk_singlets <- pred_results[pred_results$first_type == '4' & pred_results$spot_class == 'singlet',]
purk_weights <- RCTDpred@results$weights_doublet[rownames(purk_singlets),'first_type']
names(purk_weights) <- purk_singlets$singlet_score
purk_weights <- as.data.frame(purk_weights)
purk_weights$singlet_score <- as.numeric(rownames(purk_weights))

my_pal = pals::coolwarm(20)
coeff <- 1000
p2 <- ggplot(purk_weights, aes(x=singlet_score)) + 
  geom_histogram(bins=15, fill=my_pal[2], col='#ffffff', alpha=0.9) +
  stat_summary_bin(aes(singlet_score, (purk_weights-0.8)*coeff),fun='mean', bins = 14, size=0.5, geom='line') +
  stat_summary_bin(aes(singlet_score, (purk_weights-0.8)*coeff),fun='mean', bins = 14, size=1, geom='point') +
  scale_y_continuous(
    name = "Count",
    sec.axis = sec_axis(~./coeff + 0.8, name="Purkinje Proportion")
  ) +           
  theme_classic()
p2


# singlet_score - min_score vs. purkinje proportion
gen_singlet_min_plot <- function(clusters, coeff) {
  purk_singlets <- pred_results[pred_results$first_type %in% clusters & pred_results$spot_class == 'singlet',]
  purk_weights <- RCTDpred@results$weights_doublet[rownames(purk_singlets),'first_type']
  names(purk_weights) <- purk_singlets$singlet_score - purk_singlets$min_score
  purk_weights <- as.data.frame(purk_weights)
  purk_weights$singlet_score <- as.numeric(rownames(purk_weights))
  
  my_pal = pals::coolwarm(20)
  coeff <- coeff
  p2 <- ggplot(purk_weights, aes(x=singlet_score)) + 
    geom_histogram(bins=26, fill=my_pal[2], col='#ffffff', alpha=0.9) +
    stat_summary_bin(aes(singlet_score, (purk_weights-0.5)*coeff),fun='mean', bins = 26, size=0.5, geom='line') +
    stat_summary_bin(aes(singlet_score, (purk_weights-0.5)*coeff),fun='mean', bins = 26, size=1, geom='point') +
    scale_y_continuous(
      name = "Count",
      sec.axis = sec_axis(~./coeff + 0.5, name="Proportion")
    ) +           
    theme_classic() +
    xlab('singlet_score - min_score')
  return(p2)
}


gen_singlet_min_plot(c('6'), 400)


## singlet score vs marker gene
# singlet score vs. marker gene expression
marker_data_de <- readRDS(file.path(datadir, '/Data/marker_data_de_standard.RDS'))
berg_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Bergmann"]
purk_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Purkinje"]
purk_singlets <- pred_results[pred_results$first_type == '4' & pred_results$spot_class == 'singlet',]

purk_puck <- puck@counts[purk_genes,rownames(purk_singlets)]
counts_df <- as.data.frame(colSums(purk_puck))
counts_df$singlet_score <- purk_singlets$singlet_score
colnames(counts_df) <- c('counts', 'score')
p2 <- ggplot(counts_df, aes(x=score, y=counts)) + geom_point(size = 0.5)
p2

p2 <- ggplot(counts_df, aes(x=score, y=counts)) + 
  stat_summary_bin(aes(score, counts),fun='mean', bins = 8, size=0.3, geom='line') +
  stat_summary_bin(aes(score, counts),fun='mean', bins = 8, size=0.8, geom='point') +
  stat_summary_bin(aes(score, counts),fun.data = mean_se, geom='errorbar',
                   bins= 8, size= 0.3)  +
  theme_classic()
p2


# singlet_score - min_score vs. marker gene expression
marker_data_de <- readRDS(file.path(datadir, '/Data/marker_data_de_standard.RDS'))
berg_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Bergmann"]
purk_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Purkinje"]
purk_singlets <- pred_results[pred_results$first_type == '4' & pred_results$spot_class == 'singlet',]

purk_puck <- puck@counts[purk_genes,rownames(purk_singlets)]
counts_df <- as.data.frame(colSums(purk_puck))
counts_df$singlet_score <- purk_singlets$singlet_score - purk_singlets$min_score
colnames(counts_df) <- c('counts', 'score')
p2 <- ggplot(counts_df, aes(x=score, y=counts)) + geom_point(size = 0.5)
p2

p2 <- ggplot(counts_df, aes(x=score, y=counts)) + 
  stat_summary_bin(aes(score, counts),fun='mean', bins = 6, size=0.3, geom='line') +
  stat_summary_bin(aes(score, counts),fun='mean', bins = 6, size=0.8, geom='point') +
  stat_summary_bin(aes(score, counts),fun.data = mean_se, geom='line',
                   bins= 6, size= 0.3)  +
  theme_classic() +
  xlab('singlet_score - min_score')
p2


# VARIABLE UMI ANALYSIS
# doublet proportion accuracy for different total UMIs
df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df) <-  c('true', 'means', 'stdev', 'nUMI')

nUMI_vec = c(100, 200, 300, 400, 500)
for (UMI in nUMI_vec) {
  fixedUMI_results <- results[results$nUMI == UMI,]
  fixedUMI_pred_results <- pred_results[rownames(fixedUMI_results),]
  
  berg_purk_truth <- fixedUMI_results[fixedUMI_results$type1 == 'Astrocytes' & fixedUMI_results$type2 == 'Oligodendrocytes',]
  berg_purk_pred <- fixedUMI_pred_results[rownames(berg_purk_truth),]
  
  #berg_purk_pred_doublets <- berg_purk_pred[(berg_purk_pred$spot_class == 'doublet_certain') &
  #                                            ((berg_purk_pred$first_type == '0' & berg_purk_pred$second_type == '7') | (berg_purk_pred$first_type == '7' & berg_purk_pred$second_type == '0')),]
  #berg_first_IDs <- rownames(berg_purk_pred_doublets[berg_purk_pred_doublets$first_type == '7',])
  #berg_second_IDs <- rownames(berg_purk_pred_doublets[berg_purk_pred_doublets$second_type == '7',])
  
  berg_first_IDs <- rownames(berg_purk_pred[berg_purk_pred$first_type == '7',])
  berg_second_IDs <- rownames(berg_purk_pred[berg_purk_pred$second_type == '7',])
  
  berg_first_truth <- strtoi(berg_purk_truth[berg_first_IDs, 'UMI1']) / UMI
  berg_second_truth <- strtoi(berg_purk_truth[berg_second_IDs, 'UMI1']) / UMI
  berg_first_pred <- RCTDpred@results$weights_doublet[berg_first_IDs, 'first_type']
  berg_second_pred <- RCTDpred@results$weights_doublet[berg_second_IDs, 'second_type']
  berg_true_prop <- c(berg_first_truth, berg_second_truth)
  berg_pred_prop <- c(berg_first_pred, berg_second_pred)
  prop_df <- data.frame(factor(berg_true_prop), berg_pred_prop)
  true_proportions <- levels(prop_df$factor.berg_true_prop.)
  means <- c()
  stdevs <- c()
  for (prop in true_proportions) {
    means <- append(means, mean(prop_df[prop_df$factor.berg_true_prop. == prop,]$berg_pred_prop))
    stdevs <- append(stdevs, sd(prop_df[prop_df$factor.berg_true_prop. == prop,]$berg_pred_prop))
  }
  results_df <- data.frame(as.numeric(true_proportions), means, stdevs)
  colnames(results_df) <- c('true', 'means', 'stdev')
  results_df['nUMI'] <- as.character(UMI)
  df <- rbind(df, results_df)
}


my_pal = pals::coolwarm(20)
p1 <-ggplot2::ggplot(df, ggplot2::aes(x=true, y=means, group=nUMI, col=nUMI)) +
  ggplot2::geom_line() +
  ggplot2::geom_point()+
  ggplot2::geom_line(ggplot2::aes(y=true,colour = "Truth")) +
  #ggplot2::geom_errorbar(ggplot2::aes(ymin=means-1.96*stdevs, ymax=means+1.96*stdevs), width=.05) + 
  theme_classic() + xlab('True Astrocyte Proportion')+ ylab('Predicted Astrocyte Proportion') + 
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-.03,1.03)) + 
  scale_x_continuous(breaks = c(0,0.5,1), limits = c(-.03,1.03))
p1




# doublet classification rate for different total UMIs
df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df) <-  c('prop', 'value', 'std', 'nUMI')

results$UMI1 <- as.numeric(results$UMI1); results$UMI2 <- as.numeric(results$UMI2)
results$nUMI <- results$UMI1 + results$UMI2
nUMI_vec = c(100, 200, 300, 400, 500)
for (UMI in nUMI_vec) {
  fixedUMI_results <- results[results$nUMI == UMI,]
  fixedUMI_pred_results <- pred_results[rownames(fixedUMI_results),]
  doublets_pred <- fixedUMI_pred_results[fixedUMI_pred_results$spot_class == 'doublet_certain',]
  doublets_true <- fixedUMI_results[rownames(doublets_pred),]
  fixedUMI_results$prop <- factor(pmin(fixedUMI_results$UMI1, fixedUMI_results$UMI2) / UMI)
  doublets_true$prop <- factor(pmin(doublets_true$UMI1, doublets_true$UMI2) / UMI, levels=levels(fixedUMI_results$prop))
  doublet_classification_rate <- table(doublets_true$prop) / table(fixedUMI_results$prop)
  
  doublet_df <- as.data.frame(doublet_classification_rate)
  doublet_df$std <- sqrt(doublet_classification_rate * (1 - doublet_classification_rate) / table(fixedUMI_results$prop))
  colnames(doublet_df) <- c('prop', 'value', 'std')
  doublet_df$prop <- as.numeric(levels(doublet_df$prop))
  doublet_df$nUMI <- as.character(UMI)
  df <- rbind(df, doublet_df)
}


p2 <- ggplot(df, aes(x=prop,y=value,group=nUMI,col=nUMI)) + geom_line()  + 
  geom_errorbar(aes(ymin=value-1.96*std, ymax=value+1.96*std), width=.02) + theme_classic() + 
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + 
  geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.5) + 
  xlab("UMI Proportion of Minority Cell Type") + ylab("Doublet Classification Rate") + 
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1))+ scale_x_continuous(breaks = c(0,0.25,0.5), limits = c(-.03,0.53))
p2




# singlet_score - min_score vs. proportion for different total UMIs

#individual plots
gen_singlet_min_plot <- function(clusters, UMI, coeff) {
  fixedUMI_results <- results[results$nUMI == UMI,]
  fixedUMI_pred_results <- pred_results[rownames(fixedUMI_results),]
  purk_singlets <- fixedUMI_pred_results[fixedUMI_pred_results$first_type %in% clusters & fixedUMI_pred_results$spot_class == 'singlet',]
  purk_weights <- RCTDpred@results$weights_doublet[rownames(purk_singlets),'first_type']
  names(purk_weights) <- purk_singlets$singlet_score - purk_singlets$min_score
  purk_weights <- as.data.frame(purk_weights)
  purk_weights$singlet_score <- as.numeric(rownames(purk_weights))
  
  my_pal = pals::coolwarm(20)
  coeff <- coeff
  p2 <- ggplot(purk_weights, aes(x=singlet_score)) + 
    geom_histogram(bins=26, fill=my_pal[2], col='#ffffff', alpha=0.9) +
    stat_summary_bin(aes(singlet_score, (purk_weights-0.5)*coeff),fun='mean', bins = 26, size=0.5, geom='line') +
    stat_summary_bin(aes(singlet_score, (purk_weights-0.5)*coeff),fun='mean', bins = 26, size=1, geom='point') +
    scale_y_continuous(
      name = "Count",
      sec.axis = sec_axis(~./coeff + 0.5, name="Proportion")
    ) +           
    theme_classic() +
    xlab('singlet_score - min_score')
  return(p2)
}

gen_singlet_min_plot(c('0','1','2','3','4','5','6','7','8'), 500, 3000)

#aggregate plots
df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <-  c('purk_weights','singlet_score','nUMI')
for (UMI in c(100,200,300,400,500)) {
  clusters = c('0','1','2','3','4','5','6','7','8')
  fixedUMI_results <- results[results$nUMI == UMI,]
  fixedUMI_pred_results <- pred_results[rownames(fixedUMI_results),]
  purk_singlets <- fixedUMI_pred_results[fixedUMI_pred_results$first_type %in% clusters & fixedUMI_pred_results$spot_class == 'singlet',]
  purk_weights <- RCTDpred@results$weights_doublet[rownames(purk_singlets),'first_type']
  names(purk_weights) <- purk_singlets$singlet_score - purk_singlets$min_score
  purk_weights <- as.data.frame(purk_weights)
  purk_weights$singlet_score <- as.numeric(rownames(purk_weights))
  purk_weights$nUMI <- as.character(UMI)
  df <- rbind(df, purk_weights)
}

my_pal = pals::coolwarm(20)
coeff <- 15000
p2 <- ggplot(df, aes(x=singlet_score, group=nUMI, col=nUMI)) + 
  geom_histogram(bins=26, fill=my_pal[2], col='#ffffff', alpha=0.9) +
  stat_summary_bin(aes(singlet_score, (purk_weights-0.5)*coeff, group=nUMI, col=nUMI),fun='mean', bins = 26, size=0.5, geom='line') +
  stat_summary_bin(aes(singlet_score, (purk_weights-0.5)*coeff, group=nUMI, col=nUMI),fun='mean', bins = 26, size=1, geom='point') +
  scale_y_continuous(
    name = "Count",
    sec.axis = sec_axis(~./coeff + 0.5, name="Proportion")
  ) +         
  geom_hline(yintercept=0, linetype="dashed", color = "grey", size=0.5) + 
  geom_hline(yintercept=(1-0.5)*coeff, linetype="dashed", color = "grey", size=0.5) + 
  theme_classic() +
  xlab('singlet_score - min_score')
p2
