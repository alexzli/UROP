datadir <- '~/UROP'
source(file.path(datadir, '/R/analysis.R'))

RCTD_list <- readRDS(file.path(datadir, 'Objects/cerebellum_5_15.rds'))
RCTD_pred <- RCTD_list[[length(RCTD_list)]]
RCTD_truth <- readRDS(file.path(datadir, 'data/RCTD_cerebellum_slideseq.rds'))
pred_results <- RCTD_pred@results$results_df
truth_results <- RCTD_truth@results$results_df

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

marker_data_de <- readRDS(file.path(datadir, 'Data/marker_data_de_standard.RDS'))
berg_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Bergmann"]
purk_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Purkinje"]

berg_log_df <- log(RCTD_pred@cell_type_info$renorm[[1]][berg_genes,c('2','5')])
purk_log_df <- log(RCTD_pred@cell_type_info$renorm[[1]][purk_genes,c('2','5')])
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
    geom_abline(intercept = 0, slope = 1, color="black", size=0.5, linetype="dashed") + xlim(c(-15, -3)) + ylim(c(-15, -3))
p1

