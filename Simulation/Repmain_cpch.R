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

Sys.setenv(RETICULATE_PYTHON="/home/ninht/virtualenv/python3.10.4/bin/python")

setwd("Python")
library(reticulate)
source_python("cpch_source.py")

cmd_args <- commandArgs(TRUE)

xcoef_index = as.numeric(cmd_args[1])
mu_index = as.numeric(cmd_args[2])
u_n_index = as.numeric(cmd_args[3])
paral_index = as.numeric(cmd_args[4])

nsims <- 10
xcoef_options <- c(0.0, 1.0, 1.25)
mu_options <- c(0.74, 0.76, 0.78, 0.80, 0.82)
u_n_options <- list(c(2,2),c(2,3),c(3,3),c(3,4),c(4,4),c(3,5),c(4,5),c(5,5))
paral_options <- 1:50
rho <- 0

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
  
  PC_p_vec <- cpch(T = qnorm(p = punif(p_mat)), m = as.integer(n),
                   r = as.integer(u), f = f_fisher)
  #theta <- 25*zeta - 16.5
  #Z <- H*matrix(sample(x = c(theta,-theta), size = m*n, replace = TRUE),m,n) + matrix(rnorm(n = m*n),m,n)
  #PC_p_vec <- cpch(T = Z, m = as.integer(n), r = as.integer(u), f = f_fisher)
  #P <- p_mat <- 2*(1 -  pnorm(q = abs(Z)))
  
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
      R_set <- which(p.adjust(p = PC_p_vec, method = "BH") <= q)
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
      lambda <- 0.5
      pi_hat <- (1 + sum(PC_p_vec > lambda))/(m*(1-lambda))
      rej.index <- max(which(sort(PC_p_vec) < (1:m)*q/(m * pi_hat)))
      if(length(rej.index) == 0){
        return(rej.index)
      }
      R_set <- which(PC_p_vec <= sort(PC_p_vec)[rej.index])
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    
    if(method == "BY"){
      R_set <- which(p.adjust(p = PC_p_vec, method = "BY") <= q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    
    if(method == "CAMT"){
      R_set <- c()
      try(R_set <- which(CAMT::camt.fdr(pvals = PC_p_vec, pi0.var = X_cov, f1.var = X_cov)$fdr <= q), silent = TRUE)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    
    if(method == "AdaPT"){
      R_set <- c()
      n.covariates <- ncol(cbind(X_cov))
      x <- as.data.frame(X_cov)
      names(x) <- paste0("x", 1:n.covariates)
      try(Rset <- as.numeric(adaptMT::adapt_glm(x = x, pvals = PC_p_vec,
                                                pi_formulas = paste0("x", 1:n.covariates),
                                                mu_formulas = paste0("x", 1:n.covariates), 
                                                verbose = list(print = FALSE, fit = FALSE, ms = FALSE), alphas = q)$rejs), silent = TRUE)
      
      #R_set <- as.numeric(res$rejs)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    
    if(method == "IHW"){
      ihwRes <- IHW::ihw(PC_p_vec ~ rowMeans(X_cov), alpha = q)
      R_set <- which(IHW::adj_pvalues(ihwRes) <= q)
      rates <- find_FDP_TDP(rej = R_set, H = H, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$FDP/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$TDP/nsims
    }
    
  }
  
}

end_time <- Sys.time()

##
if(rho == 0){
  file.name <- paste("cpch/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
} else if(rho > 0){
  file.name <- paste("cpch/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
} else if(rho < 0){
  file.name <- paste("cpch/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
}
save(FDR_list, TPR_list, file = file.name)

print(end_time - start_time)



