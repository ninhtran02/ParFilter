# AdaFilter 
AdaFilter_procedure <- function(p_mat, u, q){
  AF_rej <- which(adaFilter::adaFilter(p.matrix = p_mat, r = u, type.I.err = "FDR", alpha = q)$decision == 1)
  return(AF_rej)
}

# BH procedure
BH_procedure <- function(p_mat, u, q, method = "Stouffer"){
  
  if(method == "Stouffer"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = stouffer_fun, u = u)
  }
  if(method == "Fisher"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = fisher_fun, u = u)
  }
  if(method == "Simes"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = simes_fun, u = u)
  }
  
  R_set <- which(p.adjust(p = p_values, method = "BH") <= q)
  
  return(R_set)
}

# CoFilter
CoFilter_procedure <- function(p_mat, u, q, tau = 0.10){
  
  p_values <- apply(X = p_mat, MARGIN = 1, FUN = fisher_fun, u = u)
  
  S <- which(p_values <= tau)
  if(length(S) == 0){
    return(c())
  }
  
  adj_p_values <- p_values[S]/tau
  
  R_res <- which(p.adjust(p = adj_p_values, method = "BH") <= q)
  
  R_set <- S[R_res]
  
  return(R_set)
}

Adaptive_CoFilter_procedure <- function(p_mat, u, q, tau = 0.10){
  
  p_values <- apply(X = p_mat, MARGIN = 1, FUN = fisher_fun, u = u)
  
  S <- which(p_values <= tau)
  if(length(S) == 0){
    return(c())
  }
  
  
  adj_p_values <- p_values[S]/tau
  
  pi_0 <- (1 + sum(adj_p_values > 0.5))/((1 - 0.5)*length(S))
  
  R_res <- which(p.adjust(p = adj_p_values, method = "BH") <= q/pi_0)
  
  R_set <- S[R_res]
  
  return(R_set)
}

# Adaptive-BH procedure
Adaptive_BH_procedure <- function(p_mat, u, q, method = "Stouffer"){
  
  if(method == "Stouffer"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = stouffer_fun, u = u)
  }
  if(method == "Fisher"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = fisher_fun, u = u)
  }
  if(method == "Simes"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = simes_fun, u = u)
  }
  
  lambda <- 0.5
  m <- length(p_values)
  pi_hat <- (1 + sum(p_values > lambda))/(m*(1-lambda))
  rej.index <- max(which(sort(p_values) < (1:m)*q/(m * pi_hat)))
  if(length(rej.index) == 0){
    return(rej.index)
  }
  R_set <- which(p_values <= sort(p_values)[rej.index])
  return(R_set)
}

# CAMT procedure
CAMT_procedure <- function(p_mat, X_mat, u, q, method = "Stouffer"){
  
  if(method == "Stouffer"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = stouffer_fun, u = u)
  }
  if(method == "Fisher"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = fisher_fun, u = u)
  }
  if(method == "Simes"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = simes_fun, u = u)
  }
  
  R_set <- which(CAMT::camt.fdr(pvals = p_values, pi0.var = X_mat, f1.var = X_mat)$fdr <= q)
  
  return(R_set)
}

# AdaPT procedure
AdaPT_procedure <- function(p_mat, X_mat, u, q, method = "Stouffer"){
  
  if(method == "Stouffer"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = stouffer_fun, u = u)
  }
  if(method == "Fisher"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = fisher_fun, u = u)
  }
  if(method == "Simes"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = simes_fun, u = u)
  }
  
  n.covariates <- ncol(cbind(X_mat))
  x <- as.data.frame(X_mat)
  names(x) <- paste0("x", 1:n.covariates)
  res <- adaptMT::adapt_glm(x = x, pvals = p_values,
                            pi_formulas = paste0("x", 1:n.covariates),
                            mu_formulas = paste0("x", 1:n.covariates), 
                            verbose = list(print = FALSE, fit = FALSE, ms = FALSE), alphas = q)
  
  R_set <- as.numeric(res$rejs)
  
  return(R_set)
}

IHW_procedure <- function(p_mat, X_mat, u, q, method = "Stouffer"){
  
  if(method == "Stouffer"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = stouffer_fun, u = u)
  }
  if(method == "Fisher"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = fisher_fun, u = u)
  }
  if(method == "Simes"){
    p_values <- apply(X = p_mat, MARGIN = 1, FUN = simes_fun, u = u)
  }
  
  ihwRes <- IHW::ihw(p_values ~ X_mat, alpha = q)
  R_set <- which(IHW::adj_pvalues(ihwRes) <= q)
  return(R_set)
}

#repfdr_procedure <- function(z_mat, u, q){
#  
#  input.to.repfdr3 <- repfdr::ztobins(zmat = z_mat, n.association.status = 2)
#  pbz <- input.to.repfdr3$pdf.binned.z
#  bz <- input.to.repfdr3$binned.z.mat
#  H <- repfdr::hconfigs(n, n.association.status = 2, studies.names = NULL)
#  
#  output.rep <- repfdr::repfdr(pbz, bz, non.null = "user.defined", non.null.rows = which(rowSums(H) >= u))
#  BayesFdr.rep <- output.rep$mat[,"Fdr"]
#  Rej <- (BayesFdr.rep <= q)
#  R_set <- which(Rej)
#  return(R_set)
#}
