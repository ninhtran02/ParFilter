# COMPUTE REJECTION SET Rrep
compute_Srep <- function(p, t, u, groups, u_groups,
                         selections, method){
  # INPUTS
  # p is a matrix of p-values
  # t is a numeric of thresholds
  # groups is a list of numerics denoting which studies belong to which groups
  # u_groups is a vector of the replicability thresholds for each group
  # selections is a list of numerics denoting the selections
  # method is the combining method

  # Last modified: Ninh Tran 7 Mar 2024

  K <- length(groups) # Number of groups
  n <- dim(p)[2]      # Number of studies
  m <- dim(p)[1]      # Number of features

  S_sf <- 1:m # Rejection set so far...
  for(k in 1:K){
    p_temp <- cbind(p[,groups[[k]]])
    p_list <- split(c(t(p_temp)), rep(1:nrow(p_temp), each = ncol(p_temp)))
    u_list <- split(rep(u_groups[k],m), 1:m)
    method_list <- split(rep(method,m), 1:m)
    simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
    simes_pu <- unname(simes_pu)

    Simes_A <- intersect(which(simes_pu <= t[n+k]),selections[[k]])
    S_sf <- intersect(S_sf,Simes_A)
  }
  return(S_sf)
}

# FDP_hat_j FUNCTION
FDP_hat_j <- function(p,t,q,j,u, groups,
                      u_groups, selections,
                      adaptive, lambda, method){
  # INPUTS
  # p is a matrix of p-values
  # t is a numeric of thresholds
  # q is a numeric of the error targets
  # j is the study of interest
  # u is the replicability threshold
  # groups is a list of numerics denoting which studies belong to which groups
  # u_groups is a numeric of the replicability thresholds for each group
  # selections is a list of numerics denoting the selections
  # adaptive is a numeric indicating whether to use adaptivity or not
  # lambda is a numeric of tuning parameters for adaptivity
  # method is the combining method

  # Last modified: Ninh Tran 7 Mar 2024

  m = nrow(p)
  n = ncol(p)
  K <- length(groups)

  for(k in 1:K){
    if(j %in% groups[[k]]){
      valid_selection_k <- intersect( which(!is.na(p[,j])), selections[[k]] )
      if(adaptive[j]){
        pi_j <- (1 + sum(p[valid_selection_k,j] > lambda[j], na.rm = TRUE))/( (1 - lambda[j])*length(valid_selection_k) )
      } else {
        pi_j <- 1
      }

      numer <- length(valid_selection_k)*t[j]*pi_j
      break
    }
  }

  denom <- which(p[,j] <= t[j])
  denom <- intersect(denom, valid_selection_k)
  Srep <- compute_Srep(p = p, t = t, u = u, groups = groups,
                       u_groups = u_groups, selections = selections,
                       method = method)
  denom <- intersect(Srep, denom)
  denom <- max(1,length(denom))

  return(numer/denom)
}

FDP_hat_repk <- function(p,t,q,k,u, groups, u_groups, selections,
                         adaptive, lambda, method = method){
  # INPUTS
  # p is a matrix of p-values
  # t is a numeric of thresholds
  # q is a numeric of the error targets
  # k is the group of interest
  # u is the replicability threshold
  # groups is a list of numerics denoting which studies belong to which groups
  # u_groups is a numeric of the replicability thresholds for each group
  # selections is a list of numerics denoting the selections
  # adaptive is a numeric indicating whether to use adaptivity or not
  # lambda is a numeric of tuning parameters for adaptivity
  # method is the combining method

  # Last modified: Ninh Tran 7 Mar 2024

  m = nrow(p)
  n = ncol(p)
  K <- length(groups)

  Srep <- compute_Srep(p = p, t = t, u = u, groups = groups,
                     u_groups = u_groups, selections = selections,
                     method = method)
  denom <- max(1, length(Srep) )

  if(adaptive[n + k]){
    p_temp <- cbind(p[,groups[[k]]])
    p_list <- split(c(t(p_temp)), rep(1:nrow(p_temp), each = ncol(p_temp)))
    u_list <- split(rep(u_groups[k],m), 1:m)
    method_list <- split(rep(method,m), 1:m)
    simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
    simes_pu <- unname(simes_pu)

    valid_selection_k <- intersect(which(!is.na(simes_pu)),selections[[k]])

    pi_k <- (1 + sum(simes_pu[valid_selection_k] > lambda[n+k], na.rm = TRUE))/( (1 - lambda[n+k])*length(valid_selection_k) )
  } else {
    pi_k <- 1
  }

  numer <- length(selections[[k]])*t[n + k]*pi_k

  return(numer/denom)
}


###### Repfilter function #######
Repfilter <- function(p, error_targets, u, groups = NULL, u_groups = NULL,
                      w = NULL, selections = NULL, adaptive,
                      lambda = NULL, auto = FALSE, omega = 0.5,
                      group_options = NULL, partition = NULL, method = "Simes"){

  # Last modified: Ninh Tran 2 April 2024

  m = nrow(p)
  n = ncol(p)

  if(!auto){
    if(is.null(groups)){
      if(partition == "minimum"){
        bp <- minimum_partition(p, u)
        groups <- bp$groups
        u_groups <- bp$u_groups
        w <- bp$w
      } else if(partition == "balanced"){
        bp <- balanced_partition(p = p, u = u)
        groups <- bp$groups
        u_groups <- bp$u_groups
        w <- bp$w
      } else if(partition == "maximum"){
        bp <- maximum_partition(p, u)
        groups <- bp$groups
        u_groups <- bp$u_groups
        w <- bp$w
      }
    }
  } else {
    if(is.null(group_options)){
      group_options <- allStirlingPartions(n = n, u = u)
    }
    ap <- auto_partition(p = p, u = u,
                         error_targets = error_targets, omega = omega,
                         group_options = group_options)
    groups <- ap$groups
    u_groups <- ap$u_groups
    w <- ap$w
  }

  K <- length(groups)

  if(is.null(u_groups)){
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
  }

  if(is.null(w)){
    w <- u_groups/u
  }

  q <- c(error_targets[1:n], rep(error_targets[n + 1], K)*w)
  if(length(adaptive) == 1){
    adaptive <- rep(adaptive,n + K)
  }
  if(is.null(lambda)){
    lambda <- q
  }

  for(l in 1:(n+K)){
    if(adaptive[l]){
      lambda[l] <- q[l]
    } else {
      lambda[l] <- 1
    }
  }

  if(is.null(selections)){
    selections <- recommended_selection(p = p, q = q,
                                        u = u, groups = groups,
                                        u_groups = u_groups,
                                        adaptive = adaptive, lambda = lambda,
                                        method = method)
  } else {
    selections <- rep(list(1:m),K)
  }

  # Starting thresholds
  t <- rep(1, n + K)
  prev_t <- rep(0, n+K)

  # Keep optiming until the previous and current threshold vectors match up
  while(!all(prev_t == t ) ){
    #print(c(t))
    prev_t <- t
    # Optimize the t_j's
    for(j in 1:n){
      t[j] <- optim_j(p = p, t = t, q = q, j = j, u = u, groups = groups,
                      u_groups = u_groups, selections = selections,
                      adaptive = adaptive,
                      lambda = lambda, method = method)
    }

    for(k in 1:K){
      t[n + k] <- optim_repk(p = p, t = t, q = q, k = k,
                             u = u, groups = groups,
                             u_groups = u_groups,
                             selections = selections,
                             adaptive = adaptive,
                             lambda = lambda, method = method)
    }
  }

  Srep <- compute_Srep(p = p, t = t, u = u, groups = groups, u_groups = u_groups,
                     selections = selections, method = method)

  R_study <- matrix(FALSE, nrow = m, ncol = n)
  Studywise_Rejections <- list()
  for(j in 1:n){
    R_study[Srep,j] <- p[Srep,j] <= t[j]
    Studywise_Rejections[[j]] <- which(R_study[,j])
  }

  output <- list(Studywise_Rejections = Studywise_Rejections,
                 Replicability_Rejections = Srep,
                 q = q,
                 lambda = lambda,
                 groups = groups,
                 u_groups = u_groups,
                 w = w)
  return(output)
}

#' Par-filter
#'
#' @param p a matrix of p-values
#' @param error_targets numeric of FDR targets. The first \eqn{n} elements are the study-specific FDR. The last element is the replicability FDR.
#' @param u the replicability threshold
#' @param groups the partition of the n studies. It is a list of numerics denoting which studies belong to which groups. This can be left NULL if \code{auto = TRUE} or a type of \code{partition} is specified.
#' @param u_groups numeric of the replicability thresholds for each group. If left NULL, \code{u_groups} will automatically be generated.
#' @param w numeric of error weights. If left NULL, \code{w} will automatically be generated.
#' @param selections list of numerics denoting the selections. If left NUll, \code{selections} will be automatically be generated based on whether \code{adaptive} is \code{TRUE} or \code{FALSE}.
#' @param adaptive Boolean indicating whether to use adaptivity or not.
#' @param lambda numeric of tuning parameters for adaptivity. If left NUll, \code{lambda} will be automatically be generated based on whether \code{error_targets}.
#' @param auto Boolean indicating whether to automatically compute the study partition
#' @param omega the overfitting tuning parameter for the automatic partioning algorithm
#' @param group_options list of candidate partitions(a list) for the automatic partitioning algorithm
#' @param partition the type of partition if groups is NULL. If \code{partition = "minimum"}, a neutral partition will be created with \eqn{K = 2} groups. If \code{partition = "minimum"}, a neutral partition will be created with \eqn{K = u} groups.
#' @param method the combining method
#'
#' @return A list containing the study-specific rejections, the rejections for assessing replicability,
#'  a numeric of study-specific and group-specific FDR target levels, the partition, the local replicability levels, and the local error weights.
#' @export
#'
#' @examples obj <- Repfilter(p = DMD.pvalues, error_targets = rep(0.05, 5),
#' u = 3, groups = list(c(1,3),c(2,4)),
#' group_options = c(2,1),
#' selections = NULL, u_groups = c(2,1),
#' adaptive = TRUE, lambda = NULL,
#' auto = FALSE, omega = 0.75, w = c(0.5,0.5),
#' partition = NULL, method = "Simes")
#' print(obj)
parfilter <- function(p, error_targets, u, groups = NULL, u_groups = NULL,
                      w = NULL, selections = NULL, adaptive,
                      lambda = NULL, auto = FALSE, omega = 0.5,
                      group_options = NULL, partition = NULL, method = "Simes"){

  # INPUTS
  # p is a matrix of p-values
  # q is a numeric of the error targets
  # u is the replicability threshold
  # groups is a list of numerics denoting which studies belong to which groups
  # u_groups is a numeric of the replicability thresholds for each group
  # w is a numeric of error weights
  # selections is a list of numerics denoting the selections
  # adaptive is a numeric indicating whether to use adaptivity or not
  # lambda is a numeric of tuning parameters for adaptivity
  # auto is a Boolean indicating whether to automatically compute the study partition
  # omega is the overfitting tuning parameter for the automatic partitioning algorithm
  # group_options is a list of candidate partitions (a list) for the automatic partitioning algorithm
  # partition is the type partition
  # method is the combining method

  # Last modified: Ninh Tran 2 April 2024

  return(Repfilter (p, error_targets, u, groups = NULL, u_groups = NULL,
             w = NULL, selections = NULL, adaptive,
             lambda = NULL, auto = FALSE, omega = 0.5,
             group_options = NULL, partition = NULL, method = "Simes"))

}

# OPTIMIZATION FUNCTIONS
optim_j <- function(p, t, q, j, u, groups, u_groups,
                    selections, adaptive = adaptive,
                    lambda = lambda, method = method){
  # INPUTS
  # p is a matrix of p-values
  # t is a numeric of thresholds
  # q is a numeric of the error targets
  # u is the replicability threshold
  # groups is a list of numerics denoting which studies belong to which groups
  # u_groups is a numeric of the replicability thresholds for each group
  # selections is a list of numerics denoting the selections
  # adaptive is a numeric indicating whether to use adaptivity or not
  # lambda is a numeric of tuning parameters for adaptivity
  # method is the combining method

  # Last modified: Ninh Tran 7 Mar 2024

  m = nrow(p)
  n = ncol(p)
  K <- length(groups)

  for(k in 1:K){
    if(j %in% groups[[k]]){
      valid_selection_k <- intersect( which(!is.na(p[,j])), selections[[k]] )
      if(adaptive[j]){
        pi_j <- (1 + sum(p[valid_selection_k,j] > lambda[j], na.rm = TRUE))/( (1 - lambda[j])*length(valid_selection_k) )
      } else {
        pi_j <- 1
      }

      p.sorted <- c(0,sort(p[valid_selection_k,j]))
      numer_vec <- length(valid_selection_k)*p.sorted*pi_j
      break
    }
  }

  cutoff_threshold <- max(which( (numer_vec/q[j] <= 0:length(valid_selection_k) ) ))
  alt_cutoff_threshold <- max(which(p.sorted <= lambda[j]))
  cutoff_threshold <- min(cutoff_threshold, alt_cutoff_threshold)
  if(cutoff_threshold == -Inf){return(0)}

  while(TRUE){
    #print(cutoff_threshold)
    t[j] <- p.sorted[cutoff_threshold]
    FDP_est <- FDP_hat_j(p = p, t = t, q = q, j = j, u = u,
                         groups = groups, u_groups = u_groups,
                         selections = selections, adaptive = adaptive,
                         lambda = lambda, method = method)
    if(FDP_est <= q[j]) return(t[j])
    denom <- (numer_vec[cutoff_threshold]/FDP_est)
    cutoff_threshold <- max(which( (numer_vec[1:cutoff_threshold]/denom  <= q[j]) ))
  }
}

# OPTIMIZATION FUNCTIONS
optim_repk <- function(p, t, q, k, u, groups, u_groups,
                       selections, adaptive = adaptive,
                       lambda = lambda, method = method){
  # INPUTS
  # p is a matrix of p-values
  # t is a numeric of thresholds
  # q is a numeric of the error targets
  # u is the replicability threshold
  # groups is a list of numerics denoting which studies belong to which groups
  # u_groups is a numeric of the replicability thresholds for each group
  # selections is a list of numerics denoting the selections
  # adaptive is a numeric indicating whether to use adaptivity or not
  # lambda is a numeric of tuning parameters for adaptivity
  # method is the combining method

  # Last modified: Ninh Tran 7 Mar 2024

  m = nrow(p)
  n = ncol(p)
  K <- length(groups)

  p_temp <- cbind(p[,groups[[k]]])
  p_list <- split(c(t(p_temp)), rep(1:nrow(p_temp), each = ncol(p_temp)))
  u_list <- split(rep(u_groups[k],m), 1:m)
  method_list <- split(rep(method,m), 1:m)
  simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
  simes_pu <- unname(simes_pu)
  valid_selection_k <- intersect(which(!is.na(simes_pu)),selections[[k]])
  simes_pu <- simes_pu[valid_selection_k]
  p.sorted <- c(0,sort(simes_pu))

  if(adaptive[n + k]){
    pi_k <- (1 + sum(simes_pu > lambda[n+k], na.rm = TRUE))/( (1 - lambda[n+k])*length(valid_selection_k) )
  } else {
    pi_k <- 1
  }


  numer_vec <- length(valid_selection_k)*p.sorted*pi_k

  cutoff_threshold <- max(which( (numer_vec/q[n + k] <= 0:length(valid_selection_k) ) ))
  alt_cutoff_threshold <- max(which(p.sorted <= lambda[n + K]))
  cutoff_threshold <- min(cutoff_threshold, alt_cutoff_threshold)
  if(cutoff_threshold == -Inf){return(0)}

  while(TRUE){
    t[n + k] <- p.sorted[cutoff_threshold]
    FDP_est <- FDP_hat_repk(p = p, t = t, q = q, k = k, u = u,
                         groups = groups, u_groups = u_groups,
                         adaptive = adaptive, lambda = lambda,
                         selections = selections, method = method)
    if(FDP_est <= q[n + k]) return(t[n + k])
    denom <- (numer_vec[cutoff_threshold]/FDP_est)
    cutoff_threshold <- max(which( (numer_vec[1:cutoff_threshold]/denom  <= q[n + k]) ))
  }
}






