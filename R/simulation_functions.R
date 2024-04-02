# Combining function
combine_p <- function(p, method){
  # INPUTS
  # p is a vector of p-values
  # method is the combining method
  
  # Last modified: Ninh Tran 7 Mar 2024
  
  # Remove NA p-values
  p <- na.omit(p)
  
  if(method == "Simes"){
    simes_p <- min(length(p)*sort(p)/(1:length(p))) 
    return(min(simes_p,1)) 
  }
  if(method == "Fisher"){
    return( 1 - pchisq(q = -2*sum(log(p)), df = 2*length(p)) )
  }
  if(method == "Stouffer"){
    return( 1 - pnorm(q = sum(qnorm(p = 1 - p, mean = 0, sd = 1))/sqrt(length(p)), mean = 0, sd = 1) )
  }
}

PC_u <- function(p,u, method){
  # INPUTS
  # p is a vector of p-values
  # u is the replicability threshold
  
  # Last modified: Ninh Tran 7 Mar 2024
  
  # Remove NA p-values
  p <- na.omit(p)
  pu <- sort(p)[u:length(p)]
  
  # Return NA if p-value vector has too few elements
  if(u > length(p)) return(NA)
  
  if(method == "Simes"){
    return(combine_p(p = pu, method = "Simes")) 
  }
  if(method == "Fisher"){
    return(combine_p(p = pu, method = "Fisher")) 
  }
  if(method == "Stouffer"){
    return(combine_p(p = pu, method = "Stouffer")) 
  }
}

# Compute FDP
compute_FDP <- function(rejections, true_nulls){
  V <- length(intersect(rejections,true_nulls))
  Rv1 <- max(length(rejections),1)
  return(V/Rv1)
}

# Compute TPP
compute_TPP <- function(rejections, false_nulls){
  L <- length(intersect(rejections,false_nulls))
  Sv1 <- max(length(false_nulls),1)
  return(L/Sv1)
}

# Generate data 
generate_pvalues <- function(m, n, u, mu, rho = 0){
  
  b <- ceiling(m*0.01)
  truth <- c()
  for(j in ceiling(n/2):n){
    truthblock <- matrix(0,nrow = b, ncol = n)
    truthblock[,1:j] <- 1
    truth <- rbind(truth,truthblock)
  }
  #truthC <- cbind(1,matrix(0,nrow = b,ncol = n - 1))
  #truthA <- matrix(0,nrow = b,ncol = n)
  #truthA[,1:ceiling(n/2)] <- 1
  ##truthA[,1:n] <- 1
  #truthB <- matrix(1,nrow = b,ncol = n)
  #
  #truth <- rbind(truthC,truthA,truthB)
  ##truth <- rbind(truthA,truthB)
  truth <- rbind(matrix(0,nrow = m - nrow(truth), ncol = n),truth)
  
  #truth <- matrix(0, nrow = m, ncol = n)
  #for(j in 1:n){
  #  truth[(m-b*j):m,j] <- 1
  #}
  
  rep_nonnull_indices <- which(rowSums(truth) >= u)
  study_nonnull_indices <- list()
  for(j in 1:n){
    study_nonnull_indices[[j]] <- which(truth[,j] == 1) 
  }
  
  if(rho == 0){
    Noise <- matrix(data = rnorm(n = m*n, mean = 0, sd = 1), nrow = m, ncol = n)
  }else{
    Noise <- c()
    for(j in 1:n){
       Noise <- cbind(Noise, arima.sim(list(order=c(1,0,0), ar=rho), n = m) / (1/sqrt(1 - rho^2)) )
    }
  }
  
  #p <- 1 - pnorm(truth*abs(matrix(rnorm(n = n*m, mean = mu, sd = 1), 
  #                                  nrow = m, ncol = n)) + Noise)
  p <- 1 - pnorm(truth*mu + Noise)
  p <- as.matrix(p)
  colnames(p) <- NULL
  
  rownames(p) <- NULL
  
  output <- list(rep_nonnull_indices = rep_nonnull_indices,
                 study_nonnull_indices = study_nonnull_indices,
                 p = p, truth = truth)
  return(output)
}
