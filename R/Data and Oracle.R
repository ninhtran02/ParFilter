# Two group-model

# Null probability conditional on xcoef
pi_0_fun <- function(eta0,xcoef,x_vec){
  # Last edited: 30th September 2024
  
  eta <- eta0 + xcoef*x_vec
  return(exp(eta)/(1 + exp(eta)))
}

# Create data for one study
create_study_data <- function(x_vec, eta0, xcoef,
                              effect_size, rho = 0){
  # Last edited: 30th September 2024
  
  m <- length(x_vec) # Number of features
  pi_0_vec <- pi_0_fun(eta0 = eta0,
                       xcoef = xcoef, x_vec = x_vec) # P(H_{ij} = 0| X_{ij})
  H_vec <- rbinom(n = m, size = 1, prob = 1 - pi_0_vec) # H_{ij}
  
  if(rho == 0){
    Noise <- rnorm(n = m, mean = 0, sd = 1)
  }else{
    Noise <- arima.sim(list(order=c(1,0,0), ar=rho),n = m)/(1/sqrt(1 - rho^2)) 
  }
  
  p_vec <- pnorm(Noise) #runif(n = length(H_vec))
  p_vec[which(H_vec == 1)] <- qbeta(p = p_vec[which(H_vec == 1)], shape1 = 1 - effect_size, shape2 = 7)
  
  output <- list(p_vec = p_vec, x_vec = x_vec, H_vec = H_vec)
  return(output)
}

# Create data for meta study
create_meta_study_data <- function(m, n, eta0_vec, xcoef_vec,
                                   effect_size_vec,
                                   rho_vec = 0){
  # Last edited: 30th September 2024
  
  # Make the parameters into vectors if they are not so already
  if( length(eta0_vec) == 1 ){
    eta0_vec <- rep(eta0_vec,n)
  }
  if( length(xcoef_vec) == 1 ){
    xcoef_vec <- rep(xcoef_vec,n)
  }
  if( length(effect_size_vec) == 1 ){
    effect_size_vec <- rep(effect_size_vec,n)
  }
  if( length(rho_vec) == 1 ){
    rho_vec <- rep(rho_vec, n)
  }
  
  p_mat <- c() # Matrix of p-values
  #z_mat <- c() # matrix of z-values
  x_mat <- mvtnorm::rmvnorm(n = m, rep(0,n), sigma = matrix(1.0,n,n) - diag(rep(1.0,n)) + diag(n) ) # Matrix of covariates
  H_mat <- c() # Matrix of Hypothesis indicators
  
  for(j in 1:n){
    x_vec <- x_mat[,j]
    #x_mat <- cbind(x_mat, x_vec)
    study_data <- create_study_data(x_vec = x_vec, eta0 = eta0_vec[j],
                                    xcoef = xcoef_vec[j],
                                    effect_size = effect_size_vec[j],
                                    rho = rho_vec[j])
    p_mat <- cbind(p_mat, study_data$p_vec)
    H_mat <- cbind(H_mat, study_data$H_vec)
    #z_mat <- cbind(z_mat, study_data$z_vec)
  }
  
  output <- list(p_mat = p_mat,
                 #z_mat = z_mat,
                 x_mat = x_mat,
                 H_mat = H_mat,
                 eta0_vec = eta0_vec,
                 xcoef_vec = xcoef_vec,
                 effect_size_vec = effect_size_vec,
                 rho_vec = rho_vec)
  return(output)
}

# The oracle procedure
oracle_procedure <- function(meta_study, u, alpha){
  
  m <- nrow(meta_study$p_mat) # Number of features
  n <- ncol(meta_study$p_mat) # Number of studies
  lfdr_mat <- c() # Matrix of (base) lfdr values
  for(j in 1:n){
    pi_0_vec <- pi_0_fun(eta0 = meta_study$eta0_vec[j],
                         xcoef = meta_study$xcoef_vec[j],
                         x_vec = meta_study$x_mat[,j])
    numerator <- pi_0_vec*dunif(x = meta_study$p_mat[,j]) 
    denominator <- numerator + (1-pi_0_vec)*dbeta(x = meta_study$p_mat[,j],
                                                  shape1 = 1 - meta_study$effect_size_vec[j], shape2 = 7)
    lfdr_mat <- cbind(lfdr_mat, numerator/denominator)
  }
  
  all_combinations <- expand.grid(rep(list(c(0, 1)), n))
  valid_combinations <- all_combinations[rowSums(all_combinations) <= u-1, ]
  
  ulfdr_vec <- rep(0, m) # Vector replicability lfdr values
  for( i in 1:nrow(valid_combinations) ){
    ulfdr_part_zeros <- apply(X = cbind(lfdr_mat[, which(valid_combinations[i,] == 0)  ]), 
                              MARGIN = 1, FUN = prod )
    ulfdr_part_ones <- apply(X = cbind(1- lfdr_mat[, which(valid_combinations[i,] == 1) ]),
                             MARGIN = 1, FUN = prod )
    ulfdr_vec <- ulfdr_vec + ulfdr_part_zeros*ulfdr_part_ones
  }
  
  # Determine which features to reject
  rej_index <- max(which(cumsum(sort(ulfdr_vec))/(1:m) <= alpha))
  rejections <- which(ulfdr_vec <= sort(ulfdr_vec)[rej_index] ) 
  
  return(rejections)
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


# Compute FDP and TPP
compute_fdp_tpp <- function(rejections, H, u){
  
  fd <- length(intersect(which(rowSums(H) < u),rejections))
  td <- length(intersect(which(rowSums(H) >= u),rejections))
  fdp <- fd/max(length(rejections),1)
  tpp <- td/max(length(which(rowSums(H) >= u)),1)
  
  return(list(fdp = fdp, tpp = tpp))
}