cell_type_means = RCTD@cell_type_info$renorm[[1]]
nUMI = RCTD@spatialRNA@nUMI
beads = t(as.matrix(RCTD@spatialRNA@counts[RCTD@internal_vars$gene_list_reg,]))
gene_list <- RCTD@internal_vars$gene_list_reg

S = data.matrix(cell_type_means[gene_list,]*nUMI[i])
B = beads[i,]
nUMI = nUMI[i]


solution <- numeric(dim(S)[2])
solution[] <- 1/length(solution)
names(solution) <- colnames(S)

set_global_Q_all()
set_likelihood_vars(Q_mat, X_vals)

solveWLS<-function(S,B,initialSol, nUMI, bulk_mode = F, constrain = F){
  solution<-pmax(initialSol,0)
  prediction = abs(S%*%solution)
  threshold = max(1e-4, nUMI * 1e-7)
  prediction[prediction < threshold] <- threshold
  gene_list = rownames(S)
  derivatives <- get_der_fast(S, B, gene_list, prediction, bulk_mode = bulk_mode)
  d_vec <- -derivatives$grad
  D_mat <- psd(derivatives$hess)
  norm_factor <- norm(D_mat,"2")
  D_mat <- D_mat / norm_factor
  d_vec <- d_vec / norm_factor
  epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
  A<-cbind(diag(dim(S)[2]))
  bzero<- (-solution)
  alpha = 0.3
  if(constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1 - sum(solution),bzero)
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A_const,b_const,meq=1)$solution
  } else {
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
}




RCTD <- create.object(
  restrict_puck(puck, barcodes), 
  cell_type_info, 
  gene_list_reg = gene_list_tot
)
RCTD <- run.algorithm(
  RCTD,
  doublet_mode = 'full',
  max_iter = max_iter,
  convergence_thresh = convergence_thresh
)
