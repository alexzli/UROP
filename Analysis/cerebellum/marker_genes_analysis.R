datadir <- '/UROP/'
source(file.path(datadir, '/R/analysis.R'))

RCTDtruth <- readRDS(file.path(datadir, '/Objects/RCTD_cerebellum_slideseq.rds'))
RCTDlist <- readRDS('../UROP/results/sim_puck_results_mean_init.rds')
RCTDpred <- RCTD_list[[length(RCTD_list)]]
cell.table(RCTDtruth, RCTDpred)

marker_data_de <- readRDS(file.path(datadir, '/Data/marker_data_de_standard.RDS'))
berg_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Bergmann"]
purk_genes <- rownames(marker_data_de)[marker_data_de$cell_type == "Purkinje"]

berg_log_df <- log(RCTDpred@cell_type_info$renorm[[1]][berg_genes,c('8','4')])
purk_log_df <- log(RCTDpred@cell_type_info$renorm[[1]][purk_genes,c('8','4')])
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

