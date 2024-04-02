# COMPUTE THE LOSS FUNCTION
loss_function <- function(p, u, groups, error_targets, omega){
  # INPUTS
  # p is a matrix of p-values
  # groups is a list of numerics denoting which studies belong to which groups
  # q is a vector of study-specific FDR targets
  # omega is a tuning parameter between 0 and 1
  
  # Last modified: Ninh Tran 2 April 2024
  
  m = nrow(p)
  n = ncol(p)
  K <- length(groups)
  u_groups <- rep(1,K)
  l <- 1
  while(sum(u_groups) != u){
    #print(u_groups)
    if(l > K){l <- 1}
    if(length(groups[[l]]) < u_groups[l] + 1){
      l <- l + 1
      next
    }
    u_groups[l] <- u_groups[l] + 1
    #print(u_groups)
    l <- l + 1
  }
  w <- u_groups/u
  q <- c(error_targets[1:n], rep(error_targets[n + 1], K)*w)
  
  quantityA <- compute_Srep(p = p, t = q, u = u, groups = groups, u_groups = u_groups, 
               selections = rep(list(1:m),K), method = "Stouffer")
  return(length(quantityA))
  
}

# AUTO-PARTITIONING ALGORITHM
auto_partition <- function(p, u, error_targets, omega, group_options){
  # INPUTS
  # p is a matrix of p-values
  # u is the replicability threshold
  # error_targets is a vector of study-specific FDR and FDRrep targets
  # omega is a tuning parameter between 0 and 1
  # groups_options is a list of lists of group partitions
  
  # Last modified: Ninh Tran 2 April 2024
  
  m <- dim(p)[1]
  n <- dim(p)[2]
  
  q <- error_targets
  
  tracker <- c()
  for(groups in group_options){
    if(!is.list(groups)){
      groups <- list(groups)
    }
    tracker <- c(tracker, loss_function(p = p, u = u, groups = groups,
                                        error_targets = error_targets, omega = omega))
    #print(tracker)
  }
  
  pick <- which.min(abs(tracker - quantile(x = tracker, probs = omega)))
  
  groups <- group_options[[pick]]
  if(!is.list(groups)){
    groups <- list(groups)
  }
  K <- length(groups)
  
  u_groups <- rep(1,K)
  l <- 1
  while(sum(u_groups) != u){
    #print(u_groups)
    if(l > K){l <- 1}
    if(length(groups[[l]]) < u_groups[l] + 1){
      l <- l + 1
      next
    }
    u_groups[l] <- u_groups[l] + 1
    #print(u_groups)
    l <- l + 1
  }
  
  w <- u_groups/u
  #w <- (u_groups^2)/sum(u_groups^2)
  
  output <- list(groups = groups, u_groups = u_groups, w = w)
  return(output)
}

#########################################################################
#################### FIXED APPROACHES ###################################
#########################################################################
# MINIMUM PARTITIONING
minimum_partition <- function(p, u){
  # INPUTS
  # p is a matrix of p-values
  # u is the replicability threshold
  
  # Last modified: Ninh Tran 2 April 2024
  
  m = nrow(p)
  n = ncol(p)
  
  K <- 2
  #group_size <- rep(1,K)*floor(n/K)
  #l <- 1
  #while(sum(group_size) != n){
  #  if(l > K){l <- 1}
  #  group_size[l] <- group_size[l] + 1
  #}
  #groups <- list(1:ceiling(n/2),(ceiling(n/2)+1):n)
  groups <- split(1:n, 1:2)
  u_groups <- rep(1,K)
  l <- 1
  while(sum(u_groups) != u){
    #print(u_groups)
    if(l > K){l <- 1}
    if(length(groups[[l]]) < u_groups[l] + 1){
      l <- l + 1
      next
    }
    u_groups[l] <- u_groups[l] + 1
    #print(u_groups)
    l <- l + 1
  }
  
  w <- u_groups/u
  
  output <- list(groups = groups, u_groups = u_groups, w = w)
  return(output)
}

# BALANCED PARTITIONING
balanced_partition <- function(p, u){
  # INPUTS
  # p is a matrix of p-values
  # u is the replicability threshold
  
  # Last modified: Ninh Tran 2 April 2024
  
  m = nrow(p)
  n = ncol(p)
  
  if(u == 2){
    return(minimum_partition(p, u))
  }
  
  K <- floor(u/2)
  
  group_size <- rep(1,K)
  k <- 1
  while(sum(group_size) != n){
    if(k > K){k <- 1}
    group_size[k] <- group_size[k] + 1
    k <- k + 1
  }
  
  groups <- list()
  l <- 0
  for(k in 1:K){
    members <- c()
    for(h in 1:group_size[k]){
      l <- l + 1
      members <- c(members, l)
    }
    groups <- c(groups, list(members))
  }
  
  u_groups <- rep(1,K)
  l <- 1
  while(sum(u_groups) != u){
    #print(u_groups)
    if(l > K){l <- 1}
    if(length(groups[[l]]) < u_groups[l] + 1){
      l <- l + 1
      next
    }
    u_groups[l] <- u_groups[l] + 1
    #print(u_groups)
    l <- l + 1
  }
  
  w <- u_groups/u
  
  output <- list(groups = groups, u_groups = u_groups, w = w)
  return(output)
}

# MAXIMUM PARTITIONING
maximum_partition <- function(p, u){
  # INPUTS
  # p is a matrix of p-values
  # u is the replicability threshold
  
  # Last modified: Ninh Tran 2 April 2024
  
  m = nrow(p)
  n = ncol(p)
  K <- u
  
  groups <- split(1:n, (1:n) %% u)
  groups <- rev(groups[order(sapply(groups, length))])
  
  u_groups <- rep(1,K)
  
  w <- u_groups/u
  
  output <- list(groups = groups, u_groups = u_groups, w = w)
  
  return(output)
}

# Recommended SELECTIONS
recommended_selection <- function(p, q, u, groups, u_groups, 
                               adaptive, lambda, method){
  # INPUTS
  # p is a matrix of p-values
  # u is the replicability threshold
  # error targets is a n+1 numeric of error targets
  # method is the combining method
  
  # Last modified: Ninh Tran 2 April 2024
  
  m = nrow(p)
  n = ncol(p)
  K <- length(groups)
  
  selections <- list()
  for(k in 1:K){
    S_sf <- 1:m
    for(l in (1:K)[-k]){
      p_temp <- cbind(p[,groups[[l]]])
      p_list <- split(c(t(p_temp)), rep(1:nrow(p_temp), each = ncol(p_temp)))
      u_local <- u_groups[l]
      u_list <- split(rep( u_local ,m), 1:m)
      method_list <- split(rep(method,m), 1:m)
      simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
      simes_pu <- unname(simes_pu)
      if(adaptive[n + l]){
        S_sf <- intersect(S_sf, which(simes_pu <= lambda[n + l]))
      } else {
        S_sf <- intersect(S_sf, which(simes_pu <= q[n + l]))
      }
    }
    selections <- c(selections, list(S_sf))
  }
  return(selections)
}


