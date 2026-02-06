# Set the library path
.libPaths("~/R/lib")

## Save personal library path as a variable
lib = .libPaths()[1]

#library("partitions", lib.loc = "~/R/lib")
#library("adaFilter")
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
u_n_options <- list(c(2,2),c(2,3),c(3,3),c(3,4),c(4,4),c(3,5),c(4,5),c(5,5))
paral_options <- 1:50
rho <- -0.8

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


methods <- c("Non-adaptive-ParFilter", "Adaptive-BH",
             "Inflated-ParFilter", "BY", "ParFilter", "BH", "Inflated-AdaFilter-BH",
             "AdaFilter-BH", "CAMT", "AdaPT", "IHW", "CoFilter-BH", "Adaptive-CoFilter-BH",
             "Naive-ParFilter")

#if(rho != 0){
#  methods <- c("ParFilter", "BH", "Inflated-AdaFilter-BH",
#               "AdaFilter-BH", "CAMT", "AdaPT", "IHW",
#               "Non-adaptive-ParFilter", "Adaptive-BH",
#               "Inflated-ParFilter", "BY", "CoFilter-BH", "Adaptive-CoFilter-BH",
#               "Naive-ParFilter")
#}

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
    if(method == "ParFilter"){
      if(u == n){
        K <- n
      }
      if(u < n){
        K <- 2
      }
      R_set <- ParFilter(P = P, u = u, X_list = X_list,
                       K = K, direction = "negative",
                       q = q, lambda_vec = rep(0.5,K))
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "Non-adaptive-ParFilter"){
      if(u == n){
        K <- n
      }
      if(u < n){
        K <- 2
      }
      R_set <- ParFilter(P = P, u = u, X_list = X_list,
                         K = K, direction = "negative",
                         q = q, lambda_vec = rep(1,K))
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "Naive-ParFilter"){
      if(u == n){
        K <- n
      }
      if(u < n){
        K <- 2
      }
      X_list_no_covar <- rep(list(rep(0,m)),n)
      R_set <- ParFilter(P = P, u = u, X_list = X_list_no_covar,
                         K = K, direction = "negative",
                         q = q, lambda_vec =  rep(0.50,K))
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "Inflated-ParFilter"){
      if(u == n){
        K <- n
      }
      if(u < n){
        K <- 2
      }
      R_set <- ParFilter(P = P, u = u, X_list = X_list,
                         K = K, direction = "negative",
                         q = q, lambda_vec = rep(1,K), inflate = TRUE)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "Inflated-AdaFilter-BH"){
      R_set <- AdaFilter_procedure(p_mat = p_mat, u = u, q = q/sum(1/(1:m)))

      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "AdaFilter-BH"){
      R_set <- AdaFilter_procedure(p_mat = p_mat, u = u, q = q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "BH"){
      R_set <- BH_procedure(p_mat = p_mat, u = u, q = q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "CoFilter-BH"){
      R_set <- CoFilter_procedure(p_mat = p_mat, u = u, q = q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "Adaptive-CoFilter-BH"){
      R_set <- Adaptive_CoFilter_procedure(p_mat = p_mat, u = u, q = q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "Adaptive-BH"){
      R_set <- Adaptive_BH_procedure(p_mat = p_mat, u = u, q = q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "BY"){
      R_set <- BH_procedure(p_mat = p_mat, u = u, q = q/sum(1/(1:m)) )
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "CAMT"){
      R_set <- c()
      try(R_set <- CAMT_procedure(p_mat = p_mat, u = u, q = q, X_mat = X_cov ), silent = TRUE)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "AdaPT"){
      R_set <- AdaPT_procedure(p_mat = p_mat, X_mat = X_cov, u = u, q = q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }

    if(method == "IHW"){
      R_set <- IHW_procedure(p_mat = p_mat, X_mat = rowMeans(X_cov), u = u, q = q)
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



