# Set the library path
set.seed(123832649)
.libPaths("~/R/lib")

## Save personal library path as a variable

lib = .libPaths()[1]
#
### Install a package
#
#install.packages("REBayes",
#                 lib = lib,
#                 repos = "https://cran.ms.unimelb.edu.au/")
#
#install.packages("ashr",
#                 lib = lib,
#                 repos = "https://cran.ms.unimelb.edu.au/")


#library("partitions", lib.loc = "~/R/lib")
library("adaFilter")

source("Data and Oracle.R")
source("ParFilter functions.R")
#source("camt.cor.func.R")

cmd_args <- commandArgs(TRUE)

xcoef_index = as.numeric(cmd_args[1])
mu_index = as.numeric(cmd_args[2])
u_n_index = as.numeric(cmd_args[3])

nsims <- 500

xcoef_options <- c(0, 1.0, 1.5)
mu_options <- c(0.74, 0.76, 0.78, 0.80, 0.82)
u_n_options <- list(c(2,2),c(2,3),c(3,3),c(3,4),c(4,4),c(3,5),c(4,5),c(5,5))

xcoef <- xcoef_options[xcoef_index]
mu <- mu_options[mu_index]
u_n <- u_n_options[[u_n_index]]
u <- u_n[1]
n <- u_n[2]

print(c(xcoef,mu,u,n))

methods <- c("Non-adaptive-ParFilter", "Adaptive-BH",
             "Inflated-ParFilter", "BY", "ParFilter", "BH", "Inflated-AdaFilter-BH",
             "AdaFilter-BH", "CAMT", "AdaPT", "IHW", "Oracle")

rho <- -0.8

if(rho != 0){
  #methods <- c("ParFilter", "BY", "Inflated-AdaFilter-BH",
  #             "AdaFilter-BH", "CAMT", "AdaPT", "IHW", "Inflated-ParFilter")
  methods <- c("ParFilter", "BH", "Inflated-AdaFilter-BH",
               "AdaFilter-BH", "CAMT", "AdaPT", "IHW", 
               "Non-adaptive-ParFilter", "Adaptive-BH",
               "Inflated-ParFilter", "BY")
}

FDR_list <- rep(list(0), length(methods))
TPR_list <- rep(list(0), length(methods))

names(FDR_list) <- methods
names(TPR_list) <- methods

m <- 5000
q <- 0.05


eta0 <- log((0.01/n)^(-1/n) - 1) #2.0
eta0_vec <- rep(eta0,n)
xcoef_vec <- rep(xcoef,n) # 1.5 is also good # xcoef_vec rep(2.5,n) very good
effect_size_vec <- rep(mu,n)
rho_vec <- rep(rho,n)

for(iter in 1:nsims){
  print(c("Working correctly: 0.0",rho))
  print(paste("Iteration:", iter, "out of", nsims))

  print(xcoef_vec)
  meta_study <- create_meta_study_data(m = m, n = n, eta0_vec = eta0_vec,
                                       xcoef_vec = xcoef_vec,
                                       effect_size_vec = effect_size_vec,
                                       rho_vec = rho_vec)
  p_mat <- meta_study$p_mat
  X_list <- list()
  for(j in 1:n){
    X_list[[j]] <- meta_study$x_mat[,j]
  }

  for(method in methods){
    if(method == "ParFilter"){
      if(u == n){
        K <- u
      }
      if(u < n){
        K <- 2
      }
      R_set <- ParFilter_FDR(p_mat = p_mat, X_list = X_list, u = u, q = q, K = K,
                             method = "Stouffer", adaptive = TRUE, cross_weights = FALSE,
                             lambdas = rep(0.50,K), inflate = FALSE)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }
    
    if(method == "Non-adaptive-ParFilter"){
      if(u == n){
        K <- u
      }
      if(u < n){
        K <- 2
      }
      R_set <- ParFilter_FDR(p_mat = p_mat, X_list = X_list, u = u, q = q, K = K,
                             method = "Stouffer", adaptive = FALSE, cross_weights = FALSE,
                             lambdas = rep(1,K), inflate = FALSE)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "Inflated-ParFilter"){
      if(u == n){
        K <- u
      }
      if(u < n){
        K <- 2
      }
      R_set <- ParFilter_FDR(p_mat = p_mat, X_list = X_list, u = u, q = q, K = K,
                             method = "Stouffer", adaptive = FALSE, cross_weights = "TRUE",
                             lambdas = rep(1,K), inflate = TRUE)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "Inflated-AdaFilter-BH"){
      R_set <- which(adaFilter(p.matrix = meta_study$p_mat, r = u,
                               alpha = q/sum(1/(1:m)))$decision == 1)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "AdaFilter-BH"){
      R_set <- which(adaFilter(p.matrix = meta_study$p_mat, r = u, alpha = q)$decision == 1)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "Oracle"){
      R_set <- oracle_procedure(meta_study = meta_study, u = u, alpha = q)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "BH"){
      R_set <- BH_procedure(p_mat = p_mat, u = u, q = q)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }
    
    if(method == "Adaptive-BH"){
      R_set <- Adaptive_BH_procedure(p_mat = p_mat, u = u, q = q)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }
    
    if(method == "BY"){
      R_set <- BH_procedure(p_mat = p_mat, u = u, q = q/sum(1/(1:5000)) )
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "CAMT"){
      R_set <- c()
      try(R_set <- CAMT_procedure(p_mat = p_mat, u = u, q = q, X_mat = rowMeans(meta_study$x_mat)), silent = TRUE)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "AdaPT"){
      R_set <- AdaPT_procedure(p_mat = p_mat, X_mat = rowMeans(meta_study$x_mat), u = u, q = q)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "IHW"){
      R_set <- IHW_procedure(p_mat = p_mat, X_mat = rowMeans(meta_study$x_mat), u = u, q = q)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }

    if(method == "repfdr"){
      R_set <- repfdr_procedure(z_mat = meta_study$z_mat, u = u, q = q)
      rates <- compute_fdp_tpp(rejections = R_set, H = meta_study$H_mat, u = u)
      FDR_list[[method]] <- FDR_list[[method]] + rates$fdp/nsims
      TPR_list[[method]] <- TPR_list[[method]] + rates$tpp/nsims
    }


  }

}



##
if(rho == 0){
  file.name <- paste("Independence/","mu",mu,"setup",u,n,"xcoef",xcoef,".RD",sep = "")
} else if(rho > 0){
  file.name <- paste("PositiveDependence/","mu",mu,"setup","xcoef",xcoef,u,n,".RD",sep = "")
} else if(rho < 0){
  file.name <- paste("NegativeDependence/","mu",mu,"setup","xcoef",xcoef,u,n,".RD",sep = "")
}
save(FDR_list, TPR_list, file = file.name)




