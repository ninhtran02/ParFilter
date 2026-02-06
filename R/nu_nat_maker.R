# nu_mat_maker under independence
#install.packages("poibin")
#nu_mat_maker <- function(selections, P, X_list, u_mat, groups_list){
#  #num_null <- 4
#  #num_alt <- 4
#  
#  K <- length(selections)
#  m <- nrow(u_mat)
#  n <- length(X_list)
#  u <- mean(rowSums(u_mat))
#  pixun_hat_mat <- matrix(NA, nrow = m, ncol = K)
#  for(k in 1:K){
#    Sk <- selections[[k]]
#    training_indices <- setdiff(1:m,Sk)
#    len <- length(training_indices)
#    
#    pix_hat_mat <- matrix(NA, nrow = m, ncol = n)
#    for(j in 1:n){
#      training_p_vec <- P[training_indices,j]
#      X_j <- X_list[[j]]
#      training_covar <- cbind(X_j)[training_indices,]
#      #f1.covar <- matrix(0, nrow = nrow(cbind(training_covar)), ncol = ncol(cbind(training_covar)))
#      camt.obj <- camt.fdr(pvals = training_p_vec,
#                           pi0.var = training_covar,
#                           f1.var = training_covar)
#      pix_hat_mat[,j] <- logistic(cbind(1,X_j)%*%camt.obj$pi0.coef)
#      # Shrinkage?
#      shr <- 1
#      pix_hat_mat[,j] <- (1 - shr) + shr*pix_hat_mat[,j]
#    }
#    # HERE
#    #if(!all(X_j == 0)){
#    #  pix_hat_mat <- 1 - H
#    #}
#    pixun_hat_mat[,k] <- apply(X = 1-pix_hat_mat, MARGIN = 1, FUN = poibin::ppoibin, kk = u-1)
#  }
#  
#  # transform? 
#  nu_mat <- (1-pixun_hat_mat)/pixun_hat_mat
#  #nu_mat <- (1 - pixun_hat_mat)
#  nu_mat <- nu_mat
#  
#  # Normalize the weights
#  for(k in 1:K){
#    nu_mat[selections[[k]],k] <- length(selections[[k]])*nu_mat[selections[[k]],k]/sum(nu_mat[selections[[k]],k])
#  }
#  return(nu_mat)
#}





# nu_mat_maker under dependence
#install.packages("poibin")
nu_mat_maker_dep <- function(selections, P, X_list, u_mat, groups_list){
  #num_null <- 4
  #num_alt <- 4
  
  K <- length(selections)
  m <- nrow(u_mat)
  n <- length(X_list)
  u <- mean(rowSums(u_mat))
  
  if(length(unique(groups_list)) != 1){
    stop("The groups are the not same for every feature!")
  } 
  
  pixucn_hat_mat <- matrix(NA, nrow = m, ncol = K)
  for(k in 1:K){
    #Sk <- selections[[k]]
    #training_indices <- setdiff(1:m,Sk)
    #len <- length(training_indices)
    
    studies <- setdiff(1:n, groups_list[[1]][[k]])
    ukc <- u - u_mat[,k]

    pix_hat_mat <- matrix(NA, nrow = m, ncol = length(studies))
    for(l in 1:length(studies)){
      j <- studies[l]
      training_p_vec <- P[,j]
      training_covar <- cbind(X_list[[j]])
      camt.obj <- camt.fdr(pvals = training_p_vec, pi0.var = training_covar, f1.var = training_covar)
      pix_hat_mat[,l] <- logistic(cbind(1,X_list[[j]])%*%camt.obj$pi0.coef)
    }
    
    for(i in 1:m){
      pixucn_hat_mat[i,k] <- poibin::ppoibin(kk = ukc[i]-1, pp = 1-pix_hat_mat[i,]) 
    }
  }
  
  # transform? 
  nu_mat <- (1-pixucn_hat_mat)/pixucn_hat_mat
  nu_mat <- nu_mat
  
  # Normalize the weights
  for(k in 1:K){
    nu_mat[selections[[k]],k] <- length(selections[[k]])*nu_mat[selections[[k]],k]/sum(nu_mat[selections[[k]],k])
  }
  return(nu_mat)
}
















#nu_mat_maker_nn <- function(selections, P, X_list, u_mat, num_null = 5, num_alt = 5, nn.size = 5){
#  #num_null <- 4
#  #num_alt <- 4
#  
#  K <- length(selections)
#  m <- nrow(u_mat)
#  n <- length(X_list)
#  u <- mean(rowSums(u_mat))
#  pixun_hat_mat <- matrix(NA, nrow = m, ncol = K)
#  for(k in 1:K){
#    Sk <- selections[[k]]
#    training_indices <- setdiff(1:m,Sk)
#    len <- length(training_indices)
#    training_PC_vec <- apply(X = P[training_indices,], MARGIN = 1, FUN = stouffer_fun, u = u) 
#    training_covar <- do.call(cbind,X_list)[training_indices,]
#    camt.obj <- camt.fdr(pvals = training_PC_vec, pi0.var = training_covar, f1.var = training_covar)
#    pixun_hat_mat[,k] <- logistic(cbind(1,do.call(cbind,X_list))%*%camt.obj$pi0.coef)
#  }
#  
#  # transform? 
#  nu_mat <- (1-pixun_hat_mat)/pixun_hat_mat
#  nu_mat <- nu_mat
#  
#  # Normalize the weights
#  for(k in 1:K){
#    nu_mat[selections[[k]],k] <- length(selections[[k]])*nu_mat[selections[[k]],k]/sum(nu_mat[selections[[k]],k])
#  }
#  return(nu_mat)
#}
#

# nu_mat_maker under independence
#install.packages("poibin")
nu_mat_maker <- function(selections, P, X_list, u_mat, groups_list){
  #num_null <- 4
  #num_alt <- 4
  
  K <- length(selections)
  m <- nrow(u_mat)
  n <- length(X_list)
  u <- mean(rowSums(u_mat))
  temp_weights_mat <- matrix(NA, nrow = m, ncol = K)
  for(k in 1:K){
    Sk <- selections[[k]]
    training_indices <- setdiff(1:m,Sk)
    len <- length(training_indices)
    
    pix_hat_mat <- matrix(NA, nrow = m, ncol = n)
    for(j in 1:n){
      training_p_vec <- P[training_indices,j]
      X_j <- X_list[[j]]
      training_covar <- cbind(X_j)[training_indices,]
      #f1.covar <- matrix(0, nrow = nrow(cbind(training_covar)), ncol = ncol(cbind(training_covar)))
      camt.obj <- camt.fdr(pvals = training_p_vec,
                           pi0.var = training_covar,
                           f1.var = training_covar)
      pix_hat_mat[,j] <- logistic(cbind(1,X_j)%*%camt.obj$pi0.coef)
      # Shrinkage?
      shr <- 1
      pix_hat_mat[,j] <- (1 - shr) + shr*pix_hat_mat[,j]
    }
    
    pixuGk_hat_vec <- temp_weights_mat[,k]
    for(i in 1:m){
      pis <- pix_hat_mat[i,groups_list[[i]][[k]]]
      uik <- u_mat[i,k]
      pixuGk_hat_vec[i] <- poibin::ppoibin(kk = uik-1, pp = 1-pis)
    }
    temp_weights_mat[,k] <- pixuGk_hat_vec
  }
  
  # transform? 
  nu_mat <- (1-temp_weights_mat)/temp_weights_mat
  nu_mat <- nu_mat
  
  # Normalize the weights
  for(k in 1:K){
    nu_mat[selections[[k]],k] <- length(selections[[k]])*nu_mat[selections[[k]],k]/sum(nu_mat[selections[[k]],k])
  }
  return(nu_mat)
}




###
nu_mat_maker <- function(selections, P, X_list, u_mat, groups_list){
  #num_null <- 4
  #num_alt <- 4
  
  K <- length(selections)
  m <- nrow(u_mat)
  n <- length(X_list)
  u <- mean(rowSums(u_mat))
  temp_weights_mat <- matrix(NA, nrow = m, ncol = K)
  for(k in 1:K){
    Sk <- selections[[k]]
    training_indices <- setdiff(1:m,Sk)
    len <- length(training_indices)
    
    pix_hat_mat <- matrix(NA, nrow = m, ncol = n)
    for(j in 1:n){
      training_p_vec <- P[training_indices,j]
      X_j <- X_list[[j]]
      training_covar <- cbind(X_j)[training_indices,]
      #f1.covar <- matrix(0, nrow = nrow(cbind(training_covar)), ncol = ncol(cbind(training_covar)))
      pvals.cutoff.vec <- seq(1e-15,1,0.01)
      for(pvals.cutoff in pvals.cutoff.vec){
        camt.obj  <- tryCatch({
          # Replace with the function you want to try
          camt.fdr(pvals = training_p_vec,
                   pi0.var = training_covar,
                   f1.var = training_covar, 
                   pvals.cutoff = pvals.cutoff)
        }, error = function(e) {
          # This executes only if an error occurs
          return(FALSE)
        })
        if(!is.logical(camt.obj)){
          break
        }
      }
      
      pix_hat_mat[,j] <- logistic(cbind(1,X_j)%*%camt.obj$pi0.coef)
      # Shrinkage?
      shr <- 1
      pix_hat_mat[,j] <- (1 - shr) + shr*pix_hat_mat[,j]
    }
    
    pixuGk_hat_vec <- temp_weights_mat[,k]
    for(i in 1:m){
      pis <- pix_hat_mat[i,groups_list[[i]][[k]]]
      uik <- u_mat[i,k]
      pixuGk_hat_vec[i] <- poibin::ppoibin(kk = uik-1, pp = 1-pis)
    }
    temp_weights_mat[,k] <- pixuGk_hat_vec
  }
  
  # transform? 
  nu_mat <- (1-temp_weights_mat)/temp_weights_mat
  nu_mat <- nu_mat
  
  # Normalize the weights
  for(k in 1:K){
    nu_mat[selections[[k]],k] <- length(selections[[k]])*nu_mat[selections[[k]],k]/sum(nu_mat[selections[[k]],k])
  }
  return(nu_mat)
}



# nu_mat_maker under dependence
#install.packages("poibin")
nu_mat_maker_dep <- function(selections, P, X_list, u_mat, groups_list){
  #num_null <- 4
  #num_alt <- 4
  
  K <- length(selections)
  m <- nrow(u_mat)
  n <- length(X_list)
  u <- mean(rowSums(u_mat))
  
  if(length(unique(groups_list)) != 1){
    stop("The groups are the not same for every feature!")
  } 
  
  pixucn_hat_mat <- matrix(NA, nrow = m, ncol = K)
  for(k in 1:K){
    #Sk <- selections[[k]]
    #training_indices <- setdiff(1:m,Sk)
    #len <- length(training_indices)
    
    studies <- setdiff(1:n, groups_list[[1]][[k]])
    ukc <- u - u_mat[,k]
    
    pix_hat_mat <- matrix(NA, nrow = m, ncol = length(studies))
    for(l in 1:length(studies)){
      j <- studies[l]
      training_p_vec <- P[,j]
      training_covar <- cbind(X_list[[j]])
      
      pvals.cutoff.vec <- seq(1e-15,1,0.01)
      for(pvals.cutoff in pvals.cutoff.vec){
        camt.obj  <- tryCatch({
          # Replace with the function you want to try
          camt.fdr(pvals = training_p_vec,
                   pi0.var = training_covar,
                   f1.var = training_covar, 
                   pvals.cutoff = pvals.cutoff)
        }, error = function(e) {
          # This executes only if an error occurs
          return(FALSE)
        })
        if(!is.logical(camt.obj)){
          break
        }
      }
      
      #camt.obj <- camt.fdr(pvals = training_p_vec, pi0.var = training_covar, f1.var = training_covar)
      pix_hat_mat[,l] <- logistic(cbind(1,X_list[[j]])%*%camt.obj$pi0.coef)
    }
    
    for(i in 1:m){
      pixucn_hat_mat[i,k] <- poibin::ppoibin(kk = ukc[i]-1, pp = 1-pix_hat_mat[i,]) 
    }
  }
  
  # transform? 
  nu_mat <- (1-pixucn_hat_mat)/pixucn_hat_mat
  nu_mat <- nu_mat
  
  # Normalize the weights
  for(k in 1:K){
    nu_mat[selections[[k]],k] <- length(selections[[k]])*nu_mat[selections[[k]],k]/sum(nu_mat[selections[[k]],k])
  }
  return(nu_mat)
}

