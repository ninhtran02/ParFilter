# Set the library path
.libPaths("~/R/lib")

## Save personal library path as a variable
lib = .libPaths()[1]

source("DataGen.R")
source("ParFilter_fun.R")
source("nu_nat_maker.R")
source("Competing Methods.R")

cmd_args <- commandArgs(TRUE)

xcoef_index = as.numeric(cmd_args[1])
mu_index = as.numeric(cmd_args[2])
u_n_index = as.numeric(cmd_args[3])
paral_index = as.numeric(cmd_args[4])

nsims <- 10
xcoef_options <- c(0.0, 1.0, 1.5)
mu_options <- c(0.74, 0.76, 0.78, 0.80, 0.82)
u_n_options <- list(c(5,5),c(6,6),c(7,7),c(8,8))
paral_options <- 1:50
rho <- 0.0

m <- 5000

paral <- paral_options[paral_index]
set.seed(123832649 + 10*round(paral))

xcoef <- gamma1 <- xcoef_options[xcoef_index]
mu <- zeta <- mu_options[mu_index]
u_n <- u_n_options[[u_n_index]]
u <- u_n[1]
n <- u_n[2]

q <- 0.05

gamma0 <- log((0.01/n)^(-1/n)-1) #1.5

print(c(gamma0,gamma1,u,n))


methods <- c("ParFilter_0", "ParFilter_1",
             "ParFilter_2", "ParFilter_3",
             "No_Covar_ParFilter_0", "No_Covar_ParFilter_1",
             "No_Covar_ParFilter_2", "No_Covar_ParFilter_3")

FDR_list <- rep(list(0), length(methods))
TPR_list <- rep(list(0), length(methods))

names(FDR_list) <- methods
names(TPR_list) <- methods

start_time <- Sys.time()
for(iter in 1:nsims){
  print(c("Working correctly: 0.0",rho))
  print(paste("Iteration:", iter, "out of", nsims))
  
  data_object <- generate_data(m = m, n = n, gamma0 = gamma0, gamma1 = gamma1, zeta = zeta)
  P <- p_mat <- data_object$P
  H <- data_object$H
  X_list <- data_object$X_list
  X_cov <- do.call(cbind, X_list)
  
  for(method in methods){
    if(method == "ParFilter_0"){
      K <- n - 0
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list,
                           K = K, direction = "negative",
                           q = q, lambda_vec = rep(0.5,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    if(method == "ParFilter_1"){
      K <- n - 1
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list,
                           K = K, direction = "negative",
                           q = q, lambda_vec = rep(0.5,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    if(method == "ParFilter_2"){
      K <- n - 2
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list,
                           K = K, direction = "negative",
                           q = q, lambda_vec = rep(0.5,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    if(method == "ParFilter_3"){
      K <- n - 3
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list,
                           K = K, direction = "negative",
                           q = q, lambda_vec = rep(0.5,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    
    if(method == "No_Covar_ParFilter_0"){
      K <- n - 0
      X_list_no_covar <- rep(list(rep(0,m)),n)
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list_no_covar,
                           K = K, direction = "negative",
                           q = q, lambda_vec =  rep(0.50,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    if(method == "No_Covar_ParFilter_1"){
      K <- n - 1
      X_list_no_covar <- rep(list(rep(0,m)),n)
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list_no_covar,
                           K = K, direction = "negative",
                           q = q, lambda_vec =  rep(0.50,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    if(method == "No_Covar_ParFilter_2"){
      K <- n - 2
      X_list_no_covar <- rep(list(rep(0,m)),n)
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list_no_covar,
                           K = K, direction = "negative",
                           q = q, lambda_vec =  rep(0.50,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    if(method == "No_Covar_ParFilter_3"){
      K <- n - 3
      X_list_no_covar <- rep(list(rep(0,m)),n)
      if(K <= n){
        R_set <- ParFilter(P = P, u = u, X_list = X_list_no_covar,
                           K = K, direction = "negative",
                           q = q, lambda_vec =  rep(0.50,K))
      }else{
        R_set <- c()
      }
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    
  }
  
}

end_time <- Sys.time()

##
if(rho == 0){
  file.name <- paste("Independence/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
} else if(rho > 0){
  file.name <- paste("PositiveDependence/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
} else if(rho < 0){
  file.name <- paste("NegativeDependence/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
}
save(FDR_list, TPR_list, file = file.name)

print(end_time - start_time)



