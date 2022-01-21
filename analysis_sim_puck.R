source('../UROP/R/optimization.R')
source('../UROP/R/initialization.R')
source('../UROP/R/clustering.R')
source('../UROP/R/analysis.R')
source('../spacexr/AnalysisCSIDE/helper_functions/de_simulation_helper.R')

reference <- readRDS('../UROP/objects/reference_RCTD_vec.rds')
puck <- generate_sim_puck(levels(reference@cell_types), rownames(reference@counts), reference)
