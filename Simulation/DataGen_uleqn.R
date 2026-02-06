# Generate covariates for the generate_data function
# Covariates are normally distributed
generate_covariates_uleqn <- function(m,n){
  X <- matrix(NA, nrow = m, ncol = n)
  for(j in 1:n){
    X[,j] <- base_x <- rnorm(n = m, mean = 0, sd = 1)
  }
  return(X)
}

# Generate data
generate_data_uleqn <- function(m,n,pi1pi1,zeta,rho=0){
  
  # Generate a matrix of covariates
  X <- generate_covariates_uleqn(m = m, n = n)
  # Generate the base null probabilities given X
  cutoff <- qnorm(p = 0.075*n, mean = 0, sd = 1)
  #base_null_prob <- pnorm(q = linear_pred)
  H <- 1*(X <= cutoff)
  H_pi1pi1 <- matrix(rbinom(n = m*n, size = 1, prob = pi1pi1),m,n)
  H <- H*H_pi1pi1
  
  Z <- matrix(NA,m,n)
  for(j in 1:n){
    Z[,j] <- arima.sim(model = list(ar = rho), n = m, sd = sqrt(1 - rho^2))
  }
  
  B <- qbeta(p = pnorm(Z), shape1 = 1 - zeta, shape2 = 7)
  P <- qunif(p = pnorm(Z))
  
  P[H == 1] <- B[H == 1]
  X_list <- split(X, col(X))
  data_object <- list(P = P,
                      H = H,
                      X_list = X_list)
}  

# Find the FDP and TDP
find_FDP_TDP <- function(rej, H, u){
  H_un <- which(rowSums(H) < u)
  A_un <- which(rowSums(H) >= u)
  FD <- length(intersect(rej,H_un))
  R <- max(length(rej),1)
  TD <- length(intersect(rej,A_un))
  TP <- max(length(A_un),1)
  FDP <- FD/R
  TDP <- TD/TP
  FDP_TDP_object <- list(FDP = FDP,
                         TDP = TDP)
  return(FDP_TDP_object)
}

find_group_FDP_TDP <- function(rej, H, groups = NULL, u_vec = NULL){
  K <- length(groups)
  m <- nrow(H)
  H_un <- c()
  for(i in 1:m){
    temp <- rep(NA,K)
    for(k in 1:K){
      if(sum(H[i,groups[[k]]]) >= u_vec[k]){
        temp[k] <-TRUE
      }
    }
    if(all(temp)){
      H_un <- c(H_un, i)
    }
  }
  
  A_un <- setdiff(1:m, H_un)
  FD <- length(intersect(rej,H_un))
  R <- max(length(rej),1)
  TD <- length(intersect(rej,A_un))
  TP <- max(length(A_un),1)
  FDP <- FD/R
  TDP <- TD/TP
  FDP_TDP_object <- list(FDP = FDP,
                         TDP = TDP)
  return(FDP_TDP_object)
}






