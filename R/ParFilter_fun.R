library(CAMT)
library(proxy)

# PC p-value Combining Functions
## Stouffer PC p-value
stouffer_fun <- function(p_vec, u){
  # Inputs
  # p_vec: vector of p-values
  # u: replicability
  
  # Return 0 if u = 0
  if(u == 0){
    return(0)
  }
  n <- length(p_vec)
  p_vec_new <- sort(p_vec)[u:n]
  PC_p <- 1 - pnorm(sum(qnorm(1 - p_vec_new))/sqrt(n-u+1))
  return(PC_p)
}

stouffer_group_fun <- function(p_vec, groups, u_vec){
  K <- length(u_vec)
  gw_PC_p_vec <- rep(NA,K)
  for(k in 1:K){
    gw_PC_p_vec[k] <- stouffer_fun(p_vec = p_vec[groups[[k]]], u = u_vec[k])
  }
  PC_p <- stouffer_fun(p_vec = gw_PC_p_vec, u = K)
  return(PC_p)
}

## Fisher PC p-value
fisher_fun <- function(p_vec, u){
  # Inputs
  # p_vec: vector of p-values
  # u: replicability
  
  # Return 0 if u = 0
  if(u == 0){
    return(0)
  }
  n <- length(p_vec)
  p_vec_new <- sort(p_vec)[u:n]
  PC_p <- 1 - pchisq(-2*sum(log(p_vec_new)), df = 2*length(p_vec_new))
  return(PC_p)
}

fisher_group_fun <- function(p_vec, groups, u_vec){
  K <- length(u_vec)
  gw_PC_p_vec <- rep(NA,K)
  for(k in 1:K){
    gw_PC_p_vec[k] <- fisher_fun(p_vec = p_vec[groups[[k]]], u = u_vec[k])
  }
  PC_p <- fisher_fun(p_vec = gw_PC_p_vec, u = K)
  return(PC_p)
}

## Simes PC p-value
simes_fun <- function(p_vec, u){
  # Inputs
  # p_vec: vector of p-values
  # u: replicability
  
  # Return 0 if u = 0
  if(u == 0){
    return(0)
  }
  n <- length(p_vec)
  p_vec_new <- sort(p_vec)[u:n]
  PC_p <- min(sort(p_vec_new)*length(p_vec_new)/(1:length(p_vec_new)))
  return(PC_p)
}

simes_group_fun <- function(p_vec, groups, u_vec){
  K <- length(u_vec)
  gw_PC_p_vec <- rep(NA,K)
  for(k in 1:K){
    gw_PC_p_vec[k] <- simes_fun(p_vec = p_vec[groups[[k]]], u = u_vec[k])
  }
  PC_p <- simes_fun(p_vec = gw_PC_p_vec, u = K)
  return(PC_p)
}


# Logistic function
logistic <- function(z){
  return(exp(z)/(exp(z) + 1))
}

# From groups and u mat
groups_list_and_u_mat_maker <- function(u, K = 2, X_list = NULL,
                                        direction = "negative", inflate = FALSE){
  
  n <- length(X_list)
  m <- nrow(cbind(X_list[[1]]))
  
  X_mat <- do.call(cbind,X_list)
  groups_list <- list()
  u_mat <- matrix(0, nrow = m, ncol = K)
  
  for(j in 1:n){
    X_study <- cbind(X_list[[j]])
    nc <- ncol(X_study)
    if(nc > 1){
      pca <- prcomp(X_study, center = TRUE, scale. = TRUE)
      z <- pca$x[, 1]   # m x 1 vector (PC1 scores)
      X_list[[j]] <- z
    }
  }
  
  if(inflate){
    x_row <- sapply(X_list, mean)
    for(i in 1:m){
      adj <- ifelse(direction == "negative", yes = 1, no = -1)
      order_x_row <- order(adj*x_row)
      rank_x_row <- rank(adj*x_row)
      # If K == n, make the groups nice
      if(K == n){
        order_x_row <- 1:n
        rank_x_row <- 1:n
      }
      groups_i <- list()
      for(k in 1:K){
        groups_i[[k]] <- order_x_row[seq(k, length(order_x_row), by = K)]
        u_mat[i,k] <- sum(rank_x_row[groups_i[[k]]] <= u)
      }
      groups_list[[i]] <- groups_i
    }
    output <- list(groups_list = groups_list,
                   u_mat = u_mat)
    return(output)
  }
  
  for(i in 1:m){
    x_row <- X_mat[i,]
    adj <- ifelse(direction == "negative", yes = 1, no = -1)
    order_x_row <- order(adj*x_row)
    rank_x_row <- rank(adj*x_row)
    # If K == n, make the groups nice
    if(K == n){
      order_x_row <- 1:n
      rank_x_row <- 1:n
    }
    groups_i <- list()
    for(k in 1:K){
      groups_i[[k]] <- order_x_row[seq(k, length(order_x_row), by = K)]
      u_mat[i,k] <- sum(rank_x_row[groups_i[[k]]] <= u)
    }
    groups_list[[i]] <- groups_i
  }
  output <- list(groups_list = groups_list,
                 u_mat = u_mat)
  return(output)
}


# Vector of group weights maker
w_vec_maker <- function(K){
  return( rep(1,K)/K )
}

PC_p_mat_maker <- function(P, groups_list, u_mat){
  m <- nrow(P)
  K <- ncol(u_mat)
  PC_p_mat <- matrix(NA,m,K)
  for(k in 1:K){
    for(i in 1:m){
      #print(groups_list[[i]][[k]])
      PC_p_mat[i,k] <- stouffer_fun(p_vec = P[i,groups_list[[i]][[k]]], u = u_mat[i,k])
    }
  }
  return(PC_p_mat)
}

selection_rule <- function(PC_p_mat, w_vec, q, lambda_vec){
  K <- length(w_vec)
  m <- nrow(PC_p_mat)
  selections <- list()
  for(k in 1:K){
    S <- 1:m
    for(l in (1:K)[-k] ){
      S <- intersect(S, which(PC_p_mat[,l] <= min(w_vec[k]*q,lambda_vec[k]) )  )
    }
    selections[[k]] <- S
  }
  return(selections)
}

#nu_mat_maker <- function(selections, P, X_list, u_mat){
#
#  m <- nrow(P)
#  n <- ncol(P)
#  K <- ncol(u_mat)
#  u <- sum(u_mat[1,])
#    
#  lfdr_mat_per_group <- list()
#    
#  for(k in 1:K){
#    not_Sk <- setdiff(x = 1:m, y = selections[[k]])
#    lfdr_mat <- c()
#    for(j in 1:n){
#      camt.obj <- camt.fdr(pvals = P[not_Sk,j], pi0.var = cbind(X_list[[j]])[not_Sk,],
#                           f1.var = cbind(X_list[[j]])[not_Sk,])
#      pi0 <- logistic(1.5*cbind(1,X_list[[j]])%*%camt.obj$pi0.coef)
#      k0 <- logistic(cbind(1,X_list[[j]])%*%camt.obj$k.coef)
#      
#      P_j_tilde <- (2*(1-k0) + k0*pi0) / (4 - 2*k0)
#      
#      lfdr_vec <- pi0/(pi0 + (1-pi0)*(1-k0)*P_j_tilde^(-k0))
#      lfdr_mat <- cbind(lfdr_mat, lfdr_vec)
#    }
#    lfdr_mat_per_group[[k]] <- lfdr_mat
#  }
#  
#  combos <- t(combn(n, u-1))
#  nu_mat <- matrix(0, nrow = m, ncol = K)
#  for(k in 1:K){
#    nu_vec_opp <- rep(0, m)
#    lfdr_mat <- lfdr_mat_per_group[[k]]
#    for(l in 1:nrow(combos)){
#      nu_vec_opp <- nu_vec_opp + rowProds(x = cbind(cbind(lfdr_mat)[,combos[l,]]))*rowProds(x = cbind(cbind(1 - lfdr_mat)[,-combos[l,]]) )
#    }
#    nu_mat[,k] <- 1 - nu_vec_opp
#  }
#  
#  # Transform weights?
#  nu_mat <- nu_mat
#  
#  # Normalize the weights
#  for(k in 1:K){
#    nu_mat[selections[[k]],k] <- length(selections[[k]])*nu_mat[selections[[k]],k]/sum(nu_mat[selections[[k]],k])
#  }
#
#  return(nu_mat)
#}

pi_hat_mat_maker <- function(nu_mat,PC_p_mat,lambda_vec,selections){
  K <- length(lambda_vec)
  m <- nrow(PC_p_mat)
  pi_hat_mat <- matrix(0, nrow = m, ncol = K)
  
  for(k in 1:K){
    Sk <- selections[[k]]
    for(i in 1:m){
      pi_hat_mat[i,k] <- (nu_mat[i,k] + sum(nu_mat[Sk,k]*(PC_p_mat[Sk,k] > lambda_vec[k])))/((1-lambda_vec[k])*length(Sk))
    }
  }
  return(pi_hat_mat)
}

R_fun <- function(t_vec, selections, PC_p_mat, nu_mat, lambda_vec, pi_hat_mat){
  K <- length(lambda_vec)
  m <- nrow(nu_mat)
  
  R_set <- 1:m
  for(k in 1:K){
    PC_p_threshold <- pmin(nu_mat[,k]*t_vec[k]/pi_hat_mat[,k], lambda_vec[k])
    R_k <- intersect(which(PC_p_mat[,k] <= PC_p_threshold),selections[[k]])
    R_set <- intersect(R_set,R_k)
  }
  
  return(R_set)
}

find_optimal_t_vec <- function(selections, PC_p_mat, nu_mat,
                               lambda_vec, pi_hat_mat, q, w_vec){
  K <- length(lambda_vec)
  t_vec_prev <- rep(0,K)
  t_vec_current <- rep(1,K)
  t_candidates_list <- split(rbind(PC_p_mat,0), col(rbind(PC_p_mat,0)))
  for(k in 1:K){
    t_candidates_list[[k]] <- sort(t_candidates_list[[k]], decreasing = TRUE)
  }
  while(!all(t_vec_prev == t_vec_current)){
    
    t_vec_prev <- t_vec_current
    
    for(k in 1:K){
      for(l in 1:length(t_candidates_list[[k]])){
        t <- t_candidates_list[[k]][l]
        t_vec_current[k] <- t
        R_set <- R_fun(t_vec = t_vec_current, selections = selections, PC_p_mat = PC_p_mat, nu_mat = nu_mat,
                       lambda_vec = lambda_vec, pi_hat_mat = pi_hat_mat)
        Sk <- selections[[k]]
        FDP_hat_k <- length(Sk)*t/max(length(R_set),1)
        if(FDP_hat_k <= w_vec[k]*q){
          t_candidates_list[[k]] <- t_candidates_list[[k]][t_candidates_list[[k]] <= t] 
          break
        }
      }
    }
  }
  return(t_vec_current)
}


####

ParFilter <- function(P, u, X_list, K, direction = "negative",
                      q = q, lambda_vec = rep(0.5,K),
                      inflate = FALSE){
  
  m <- nrow(P)
  n <- ncol(P)
  if(length(unique(c(ncol(P),K,u)))==1){
    groups_list <- replicate(m, as.list(1:n), simplify = FALSE)
    u_mat <- matrix(data = 1, nrow = m, ncol = K)
  } else {
    groups_list_and_u_mat <- groups_list_and_u_mat_maker(u = u, K = K, X_list = X_list, direction = "negative", inflate = inflate)
    u_mat <- groups_list_and_u_mat$u_mat
    groups_list <- groups_list_and_u_mat$groups_list
  }
  w_vec <- w_vec_maker(K = K)
  
  PC_p_mat <- PC_p_mat_maker(P = P, groups_list = groups_list, u_mat = u_mat)
  
  if(inflate){
    # if inflate == TRUE
    selections <- selection_rule(PC_p_mat = PC_p_mat, w_vec = w_vec, q = q/(sum(1/(1:m))), lambda_vec = lambda_vec)
  } else {
    # if inflate == FALSE
    selections <- selection_rule(PC_p_mat = PC_p_mat, w_vec = w_vec, q = q, lambda_vec = lambda_vec)
  }
 
  if(inflate){
    # if inflate == TRUE
    nu_mat <- nu_mat_maker_dep(selections = selections, P = P, X_list = X_list, groups_list = groups_list,
                              u_mat = u_mat)
  } else {
    # if inflate == FALSE
    nu_mat <- nu_mat_maker(selections = selections, P = P, X_list = X_list,
                              u_mat = u_mat, groups_list = groups_list)
  }
  
  if(inflate || all(lambda_vec == 1)){
    # if inflate == TRUE
    pi_hat_mat <- matrix(data = 1, nrow = m, ncol = K)
  } else {
    # if inflate == FALSE
    pi_hat_mat <- pi_hat_mat_maker(nu_mat = nu_mat, PC_p_mat = PC_p_mat, lambda_vec = lambda_vec, selections = selections)
  }
  
  if(inflate){
    # if inflate == TRUE
    for(k in 1:K){
      w_vec[k] <- w_vec[k]/(sum(1/(1:length(selections[[k]]))))
    }
  }
  
  t_vec <- find_optimal_t_vec(selections = selections,
                              PC_p_mat = PC_p_mat, 
                              nu_mat = nu_mat, 
                              lambda_vec = lambda_vec,
                              pi_hat_mat = pi_hat_mat, q = q, w_vec = w_vec)
  
  rej <- R_fun(t_vec = t_vec, selections = selections, PC_p_mat = PC_p_mat, nu_mat = nu_mat, lambda_vec = lambda_vec, pi_hat_mat = pi_hat_mat)
  return(rej)
}