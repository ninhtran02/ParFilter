# COMBINING FUNCTIONS (for testing null PC hypothesis)
# Fisher combining function
fisher_fun <- function(p_vec, u){
  # Last edited: 3rd October 2024

  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level

  n <- length(p_vec)
  p_subset <- sort(p_vec, decreasing = TRUE)[1:(n - u + 1)]
  len_subset <- length(p_subset)
  PC_p <- 1 - pchisq(q = -2*sum(log(p_subset)), df = 2*len_subset)

  # Output: One PC p-value
  return(PC_p)
}
# Stouffer combining function
stouffer_fun <- function(p_vec, u){
  # Last edited: 3rd October 2024

  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level

  n <- length(p_vec)
  p_subset <- sort(p_vec, decreasing = TRUE)[1:(n - u + 1)]
  len_subset <- length(p_subset)
  PC_p <- 1 - pnorm(q = sum(qnorm(p = 1 - p_subset))/sqrt(len_subset))

  # Output: One PC p-value
  return(PC_p)
}
# Simes combining function
simes_fun <- function(p_vec, u){
  # Last edited: 3rd October 2024

  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level

  n <- length(p_vec)
  p_subset <- sort(p_vec, decreasing = TRUE)[1:(n - u + 1)]
  PC_p <- min(p.adjust(p = p_subset, method = "BH"))

  # Output: One PC p-value
  return(PC_p)
}

# DEFAULT SELECTION RULES
default_selections <- function(p_mat, groups, u_groups_mat, method, q, lambdas){
  # Last edited: 22nd October 2024

  # Inputs
  ## p_mat: mxn matrix of p-values
  ## groups: List containing the study groupings
  ## u: Replicability level
  ## method: Method for combining p-values: "Fisher", "Stouffer", or "Simes"
  ## q: FDR target
  ## weights: Logical for whether to weight the PC p-values

  K <- length(groups) # Number of groups
  m <- nrow(p_mat) # Number of features
  n <- ncol(p_mat) # Number of studies
  u <- rowSums(u_groups_mat)[1]
  G_sizes <- unlist(lapply(X = groups, FUN = length)) # Group sizes
  #min_Gk <- min(G_sizes) # Smallest group size
  #u_bar <- max(u - n + min_Gk,1) # Local replicability levels (for selections)
  #u_bar <- pmax(4 - 5 + G_sizes,1)
  w_vec <- G_sizes/n # Local error weights

  # Weights for the features
  weights_list <- replicate(K, rep(1,m), simplify = FALSE)

  # Create a list of local (u_bar/|Gk|) PC p-values
  PC_p_list <- list()
  for(k in 1:K){
    if(u == n){
      if(method == "Fisher"){
        PC_p_list[[k]] <- apply(X = cbind(as.matrix(p_mat[,groups[[k]]])), MARGIN = 1, FUN = fisher_fun, u = G_sizes[k])
      }
      if(method == "Stouffer"){
        PC_p_list[[k]] <- apply(X = cbind(as.matrix(p_mat[,groups[[k]]])), MARGIN = 1, FUN = stouffer_fun, u = G_sizes[k])
      }
      if(method == "Simes"){
        PC_p_list[[k]] <- apply(X = cbind(as.matrix(p_mat[,groups[[k]]])), MARGIN = 1, FUN = simes_fun, u = G_sizes[k])
      }
    } else {
      sub_p_mat <- cbind(as.matrix(p_mat[,groups[[k]]]))
      sub_p_mat_list <- split(sub_p_mat, row(sub_p_mat))
      u_groups_mat_list <- as.list(u_groups_mat[,k])
      PC_p_list[[k]] <- mapply(fisher_fun, p_vec = sub_p_mat_list, u = u_groups_mat_list)
    }

  }

  # Threshold the the local PC p-values
  which_thresholded <- list()
  for(k in 1:K){
    which_thresholded[[k]] <- which(PC_p_list[[k]] <= pmin(weights_list[[k]]*w_vec[k]*q,lambdas[k]))
  }

  selections <- list()
  for(k in 1:K){
    # Intersect thresholded features outside of group k
    selections[[k]] <- Reduce(f = intersect, x = which_thresholded[-k])
  }

  # Output: List of default selections for each group
  return(selections)
}


# COMPUTE LOCAL REPLICABILITY LEVELS
u_groups_mat_maker <- function(groups, u, lfdr_mat){
  # Last edited: 3rd October 2024

  # Inputs
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## u: Replicability level

  n <- ncol(lfdr_mat) # Number of studies
  m <- nrow(lfdr_mat) # Number of features
  K <- length(groups) # Number of groups

  if(u == n){
    u_s <- rep(unlist(lapply(groups, length)),m)
    u_groups_mat <- matrix(data = u_s, nrow = m, ncol = K, byrow = TRUE)
    return(u_groups_mat)
  }

  if(u == K){
    u_groups_mat <- matrix(data = 1, nrow = m, ncol = K)
    return(u_groups_mat)
  }

  # Transform the weights using -log()
  study_weights_mat <- 1 - lfdr_mat

  # Compute the local replicability levels for each feature
  u_groups_mat <- t(apply(X = study_weights_mat, MARGIN = 1,
                          FUN = u_groups_maker, groups = groups, u = u))

  # Output: Matrix of local replicability levels
  return(u_groups_mat)
}

# COMPUTE LOCAL REPLICABILITY LEVELS BASED ON STUDY WEIGHTS
u_groups_maker <- function(study_weights, groups, u){
  # Last edited: 3rd October 2024

  ## Inputs
  ## study_weights: Vector of study weights
  ## groups: List of groups
  ## u: Replicability level

  n <- length(study_weights) # Number of studies
  K <- length(groups) # Number of groups

  study_weight_ranks <- rank(-study_weights) # Rank of the weights

  grouped_study_weight_ranks <- list()
  u_groups <- rep(1,K)
  for(k in 1:K){
    group_k <- study_weight_ranks[groups[[k]]]
    group_k <- c(Inf,group_k[-which.min(group_k)])
    grouped_study_weight_ranks[[k]] <- group_k
  }

  while(sum(u_groups) < u){
    k <- which.min(unlist(lapply(grouped_study_weight_ranks, min)))
    u_groups[k] <- u_groups[k] + 1
    group_k <- grouped_study_weight_ranks[[k]]
    grouped_study_weight_ranks[[k]] <- group_k[-which.min(group_k)]
  }

  # Output: Vector of local replicability levels
  return(u_groups)
}

# COMPUTE APPROXIMATE LFDR
compute_approx_lfdr <- function(p_mat, X_list, groups,
                                selections, cross_weights = FALSE){
  # Last edited: 23rd October 2024

  ## Input
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## selections: List of default selections for each group

  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(selections) # Number of groups

  if(cross_weights == "naive"){
    cross_weights <- FALSE
  }

  # Make a list of covariates (matrix) for each study
  for(j in 1:n){
    X_list[[j]] <- cbind(X_list[[j]])
  }

  # Create a matrix of estimated values for Pr(H_{ij} = 0 | x_{ij})
  # Get training data from outside the groups
  training_p_values_list <- list()
  training_covariates_list <- list()
  for(j in 1:n ){
    k_index <- which( unlist(lapply(X = groups, FUN = function(x){ j %in% x })) )

    if(cross_weights){
      training_indices <- 1:m
    } else {
      training_indices <- setdiff(1:m,selections[[k_index]])
    }

    training_p_values_list[[j]] <- p_mat[training_indices,j]
    training_covariates_list[[j]] <- X_list[[j]][training_indices,]
  }

  # Train the parameters
  pi0.coef_list <- list()
  k.coef_list <- list()
  for(j in 1:n){
    camt_obj <- CAMT::camt.fdr(pvals = training_p_values_list[[j]],
                               pi0.var = training_covariates_list[[j]],
                               f1.var = training_covariates_list[[j]])

    if(cross_weights){
      pi0.coef_list[[j]] <- camt_obj$pi0.coef
    } else {
      # Add in an adjustment factor when there is selection
      pi0.coef_list[[j]] <- 1.5*camt_obj$pi0.coef
    }

    k.coef_list[[j]]   <- 1.0*camt_obj$k.coef
  }


  # Create a matrix of estimated values for Pr(H_{ij} = 0 | x_{ij})
  # for features that were selected in group k (trained on the data from
  # studies not in group k)
  pi0_est_mat <- matrix(NA, nrow = m, ncol = n)
  k_est_mat <- matrix(NA, nrow = m, ncol = n )
  for(j in 1:n){
    X_design <- cbind(1, X_list[[ j ]])
    eta <- X_design%*%(pi0.coef_list[[j]])
    beta <- X_design%*%(k.coef_list[[j]])
    pi0_est_mat[,j] <- 1/(1 + exp(-eta))
    k_est_mat[,j] <- 1/(1 + exp(-beta))
  }


  if(cross_weights){
    p_mat_for_weights <- p_mat
    p_mat_for_weights <- -(k_est_mat*pi0_est_mat - 2*k_est_mat + 2)/(2*k_est_mat - 4)
  } else {
    p_mat_for_weights <- -(k_est_mat*pi0_est_mat - 2*k_est_mat + 2)/(2*k_est_mat - 4)
  }
  hp_mat <- (1 - k_est_mat)*p_mat_for_weights^(-k_est_mat)
  #hp_mat <- (1 - k_est_mat)*pmin(p_mat_for_weights^(-k_est_mat),100)
  lfdr_mat <- pi0_est_mat/( (1 - pi0_est_mat)*hp_mat + pi0_est_mat )

  #lfdr_mat <- 1 - k_est_mat
  #lfdr_mat <- pi0_est_mat

  pi0.coef <- c()
  k.coef <- c()
  for(k in 1:K){
    pi0.coef <- c(pi0.coef, pi0.coef_list[[k]][-1])
    k.coef <-   c(k.coef,   k.coef_list[[k]][-1]  )
  }

  # Output: List of weights for selected features in each group
  return(list(lfdr_mat = lfdr_mat,
              pi0_est_mat = pi0_est_mat,
              k_est_mat = k_est_mat,
              pi0.coef = pi0.coef,
              k.coef   = k.coef))
}

# COMPUTE WEIGHTS FOR FDR CONTROL
weights_maker <- function(lfdr_mat, groups, u_groups_mat, selections,
                          normalize = TRUE){
  # Last edited: 23rd October 2024

  #print("PROTOTYPE K")

  ## Input
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## u_groups_mat: Matrix of local replicability levels
  ## selections: List of refined selections for each group

  K <- length(groups)
  weights_list <- list() # List of weights for selected features for each group

  for(k in 1:K){
    est_ulfdr_vec <- rep(0, length(selections[[k]]))
    all_combinations <- expand.grid(rep(list(c(0, 1)), length(groups[[k]]) ))
    for(i in 1:length(selections[[k]])){
      valid_combinations <- cbind(all_combinations[rowSums(all_combinations) < u_groups_mat[ selections[[k]][i] ,k], ])

      for( l in 1:nrow(valid_combinations) ){
        est_ulfdr_part_zeros <- prod(lfdr_mat[ selections[[k]][i],    groups[[k]][which(valid_combinations[l,] == 0)]  ]    )
        est_ulfdr_part_ones <- prod(1- lfdr_mat[ selections[[k]][i] , groups[[k]][which(valid_combinations[l,] == 1)]  ]    )
        est_ulfdr_vec[i] <- est_ulfdr_vec[i] + est_ulfdr_part_zeros*est_ulfdr_part_ones
      }
    }

    # Scale the weights
    if(normalize){
      feature_weights <- 1 - est_ulfdr_vec
      feature_weights <- (feature_weights/sum(feature_weights))*length(selections[[k]])
      weights_list[[k]] <- feature_weights
    } else {
      weights_list[[k]] <- 1 - est_ulfdr_vec
    }
  }

  # Output: List of weights for selected features in each group
  return(weights_list)
}

# COMPUTE WEIGHTS FOR FDR CONTROL
cross_weights_maker <- function(lfdr_mat, groups, u_groups_mat, selections,
                          normalize = TRUE){
  # Last edited: 23rd October 2024

  #print("PROTOTYPE K")

  ## Input
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## u_groups_mat: Matrix of local replicability levels
  ## selections: List of refined selections for each group

  n <- ncol(cbind(lfdr_mat))

  weights_list <- list() # List of weights for selected features for each group

  for(k in 1:K){
    est_ulfdr_vec <- rep(0, length(selections[[k]]))
    all_combinations <- expand.grid(rep(list(c(0, 1)), n - length(groups[[k]]) ))
    oosg <- unlist(groups[-k])
    for(i in 1:length(selections[[k]])){
      valid_combinations <- cbind(all_combinations[rowSums(all_combinations) < u - u_groups_mat[ selections[[k]][i] ,k], ])

      for( l in 1:nrow(valid_combinations) ){
        est_ulfdr_part_zeros <- prod(lfdr_mat[ selections[[k]][i],    oosg[which(valid_combinations[l,] == 0)]  ]    )
        est_ulfdr_part_ones <- prod(1- lfdr_mat[ selections[[k]][i] , oosg[which(valid_combinations[l,] == 1)]  ]    )
        est_ulfdr_vec[i] <- est_ulfdr_vec[i] + est_ulfdr_part_zeros*est_ulfdr_part_ones
      }
    }

    # Scale the weights
    if(normalize){
      # Transform the weights
      feature_weights <- 1 - est_ulfdr_vec
      feature_weights <- (feature_weights/sum(feature_weights))*length(selections[[k]])
      weights_list[[k]] <- feature_weights
    } else {
      weights_list[[k]] <- 1 - est_ulfdr_vec
    }
  }

  # Output: List of weights for selected features in each group
  return(weights_list)
}

# COMPUTE PC P-VALUES FOR EACH GROUP
PC_p_maker <- function(p_mat, groups, u_groups_mat, method){
  # Last edited: 3rd October 2024

  ## Inputs
  ## p_mat: mxn matrix of p-values
  ## groups: List containing the study groupings
  ## u_groups_mat: Matrix of local replicability levels
  ## method: Method for combining p-values: "Fisher", "Stouffer", or "Simes"

  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(groups) # Number of groups

  PC_p_list <- list()
  for(k in 1:K){
    if(method == "Fisher"){
      #PC_p_list[[k]] <- mapply(fisher_fun, p_vec = p_mat[,groups[[k]]], u = u_groups_mat[,k])
      PC_p_vec <- rep(NA, m)
      for(i in 1:m){
        PC_p_vec[i] <- fisher_fun(p_vec = p_mat[i,groups[[k]]], u = u_groups_mat[i,k])
      }
      PC_p_list[[k]] <- PC_p_vec
    }
    if(method == "Stouffer"){
      #PC_p_list[[k]] <- mapply(stouffer_fun, p_vec = p_mat[,groups[[k]]], u = u_groups_mat[,k])
      PC_p_vec <- rep(NA, m)
      for(i in 1:m){
        PC_p_vec[i] <- stouffer_fun(p_vec = p_mat[i,groups[[k]]], u = u_groups_mat[i,k])
      }
      PC_p_list[[k]] <- PC_p_vec
    }
    if(method == "Simes"){
      #PC_p_list[[k]] <- mapply(simes_fun, p_vec = p_mat[,groups[[k]]], u = u_groups_mat[,k])
      PC_p_vec <- rep(NA, m)
      for(i in 1:m){
        PC_p_vec[i] <- simes_fun(p_vec = p_mat[i,groups[[k]]], u = u_groups_mat[i,k])
      }
      PC_p_list[[k]] <- PC_p_vec
    }
  }

  # Output: List of PC p-values for each group
  return(PC_p_list)
}

# COMPUTE (WEIGHTED) STOREY NULL PROPORTIONE ESTIMATE (2004)
pi0_hat_maker <- function(PC_p_list,
                    selections, weights_list, lambdas){
  # Last edited: 3rd October 2024

  ## Inputs:
  ## PC_p_list: List of PC p-values for each group
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters

  K <- length(selections) # Number of groups
  pi0_hat_vec <- rep(NA,K)
  for(k in 1:K){
    num_sum_index <- which(PC_p_list[[k]][selections[[k]]] > lambdas[k])
    numerator <- max(weights_list[[k]]) + sum(weights_list[[k]][num_sum_index])
    denominator <- length(selections[[k]])*(1 - lambdas[k])
    pi0_hat_vec[k] <- numerator/denominator
  }

  # Output: Vector of (post-selection) pi0 estimates
  return(pi0_hat_vec)
}

# COMPUTE ParFilter REJECTION SET
create_R <- function(PC_p_list, t_vec, selections,
                     weights_list, lambdas){
  # Last edited: 3rd October 2024

  ## Inputs
  ## PC_p_list: List of PC p-values for each group
  ## t_vec: Vector of rejection thresholds for each group
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters

  K <- length(PC_p_list) # Number of groups
  m <- length(PC_p_list[[1]]) # Number of features
  R_set <- 1:m # Initial rejection indices
  for(k in 1:K){
    adj_p <- PC_p_list[[k]][selections[[k]]]/weights_list[[k]]
    unadj_p <- PC_p_list[[k]][selections[[k]]]
    # Find which selections in group k to reject
    R_which <- intersect(which(adj_p <= t_vec[k]),which(unadj_p <= lambdas[k]))
    R_set <- intersect(R_set, selections[[k]][R_which])
  }

  # Output: Vector of rejected indices
  return(R_set)
}

# Compute a vector of FDP_hat_k values based on a
# vector of candidate values for the group k threshold
FDP_hat_k_vec <- function(PC_p_list, t_vec, k, tk_candidates,
                          pi0_hat_vec, selections,
                          weights_list, lambdas){

  # Last edited: 3rd October 2024

  ## Inputs
  ## PC_p_list: List of PC p-values for each group
  ## t_vec: Vector of rejection thresholds for each group
  ## k: Group of interest
  ## tk_candidates: Vector of values for tk
  ## pi0_hat_vec: Vector of (post-selection) pi0 estimates
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters

  R_sizes <- rep(NA, length(tk_candidates)) # Vector to store rejection sizes
  for(i in 1:length(tk_candidates) ){
    t_vec_temp <- t_vec
    t_vec_temp[k] <- tk_candidates[i] # Replace current tk with candidate
    R_set <- create_R(PC_p_list = PC_p_list, t_vec = t_vec_temp,
                      selections = selections, weights_list = weights_list,
                      lambdas = lambdas)
    R_sizes[i] <- length(R_set)
  }

  FDP_hat_vec <- length(selections[[k]])*pi0_hat_vec[k]*tk_candidates/pmax(R_sizes,1)

  # Output: Vector of FDP_hat values for group k given candidate tk thresholds
  return(FDP_hat_vec)
}

# COMPUTE GROUPS
groups_maker <- function(u,n,K){
  # Last edited: 3rd October 2024

  group_index <- cut(1:n, breaks = K, labels = FALSE)
  groups <- split(1:n, group_index)

  # Output: List of study groupings
  return(groups)
}

# COMPUTE LOCAL ERROR WEIGHTS
w_vec_maker <- function(groups){
  # Last edited: 3rd October 2024

  w_vec <- unlist(lapply(X = groups, FUN = length))/length(unlist(groups))

  # Output: Vector of local error weights
  return(w_vec)
}

# ParFilter ALGORITHM
#' ParFilter
#'
#' @param p_mat mxn matrix of p-values
#' @param X_list List of covariates (numeric of matrix) for each study
#' @param u Replicability level. If u is less than the number of studies, then the local replicability levels will be determined randomly.
#' @param q FDR target
#' @param K Number of groups
#' @param method Method for combining p-values: "Fisher", "Stouffer", or "Simes"
#' @param adaptive logical indicating whether to use adaptive weighted null proportion estimates or not
#' @param lambdas numeric of tuning parameters for the weighted null proportion estimates.
#' @param cross_weights Set as TRUE if the p-values are dependent within studies, otherwise leave it as FALSE.
#' @param inflate Use a harmonic number as the weighted null proportion estimates.
#'
#' @return The rejections for determining u/[n] replicability
#' @export
#'
#' @examples
ParFilter_FDR <- function(p_mat, X_list, u, q, K, method,
                          adaptive = TRUE, lambdas = rep(0.5,K),
                          cross_weights = FALSE, inflate = FALSE){
  # Last edited: 22nd October 2024

  ## Input
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study
  ## u: Replicability level
  ## q: FDR target
  ## K: Number of groups
  ## method: Method for combining p-values: "Fisher", "Stouffer", or "Simes"
  ## adaptive: logical indicating whether to use adaptive null proportion estimates or not
  ## cross_weights: Set as TRUE if the p-values are dependent within studies, otherwise leave it as FALSE.
  ## lambdas: numeric of tuning parameters for the null proportion estimates.

  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features

  # Compute the study groupings and local error_weights
  groups <- groups_maker(u = u, n = n, K = K)
  w_vec <- w_vec_maker(groups)
  #print(w_vec)

  # Make u_groups_mat
  if(u < n){
    weights_for_u <- matrix(data = runif(n = m*n, min = 0, max = 1), nrow = m, ncol = n)
  } else {
    weights_for_u <- matrix(data = rep(x = 0.5, m*n), nrow = m, ncol = n)
  }
  u_groups_mat <- u_groups_mat_maker(groups = groups, u = u,
                                     lfdr_mat = weights_for_u)

  # Form the selection rules
  if(!inflate){
    selections <- default_selections(p_mat = p_mat, groups = groups,
                                     u_groups_mat = u_groups_mat, method = method, q = q,
                                     lambdas = lambdas)
  }else{
    #print("IM HERE")
    selections <- default_selections(p_mat = p_mat, groups = groups,
                                     u_groups_mat = u_groups_mat, method = method, q = q/sum(1/(1:m)),
                                     lambdas = lambdas)
  }


  # Early exit if the selections are empty
  if(0 %in% lapply(selections, length)){
    return(c())
  }

  lfdr_obj <- compute_approx_lfdr(p_mat = p_mat, X_list = X_list,
                                  groups = groups,
                                  selections = selections, cross_weights = cross_weights)

  lfdr_mat <- lfdr_obj$lfdr_mat

  # Make u_groups_mat

  # Form the weights
  if(cross_weights == "naive"){
    weights_list <- list()
    for(k in 1:K){
      weights_list[[k]] <- rep(1,length(selections[[k]]))
    }
  }
  if(cross_weights == FALSE){
    weights_list <- weights_maker(lfdr_mat = lfdr_mat, groups = groups,
                                  u_groups_mat = u_groups_mat,
                                  selections = selections)
  }
  if(cross_weights == TRUE){
    weights_list <- cross_weights_maker(lfdr_mat = lfdr_mat, groups = groups,
                                  u_groups_mat = u_groups_mat,
                                  selections = selections)
  }
  if(cross_weights == "x"){
    weights_list <- rep(list(),K)
    for(k in 1:K){
      x_selected <- rowMeans( cbind(do.call(cbind, X_list[groups[[k]]])[selections[[k]],]))
      weights_temp <- sort( (1:length(selections[[k]]))^(2) )[rank(-x_selected)]
      weights_list[[k]] <- (weights_temp/sum(weights_temp))*length(selections[[k]])
    }
  }
  #print(weights_list)

  # Form the PC p-values
  PC_p_list <- PC_p_maker(p_mat = p_mat, groups = groups,
                          u_groups_mat = u_groups_mat, method = method)

  #print(lambdas)

  # Form the pi0_hats
  if(adaptive){
    #for( k in 1:K ){
    #  lambdas[k] <- 1 - min(mean(lfdr_obj$pi0_est_mat[selections[[k]],k]*weights_list[[k]]),0.95)
    #}
    #print(lambdas)

    pi0_hat_vec <- pi0_hat_maker(PC_p_list = PC_p_list, selections = selections,
                                 weights_list = weights_list, lambdas = lambdas)
  } else {
    pi0_hat_vec <- rep(1,K)
  }

  if(inflate){
    selection_lengths <- unlist(lapply(X = selections, FUN = length))
    print(selection_lengths)
    inflation_factors <- rep(NA,K)
    for(k in 1:K){
      inflation_factors[k] <- 1/sum(1/(1:selection_lengths[k]))
    }
    print(inflation_factors)
    print(w_vec)
    w_vec <- w_vec*inflation_factors
    print(w_vec)
  }

  #print("HERE IS WEIGHT LIST LENGTH:",length(weights_list))

  # Find the rejection thresholds
  t_vec <- rep( m, K )
  t_vec_prev <- rep(0,K)

  while(!all(t_vec == t_vec_prev)){
    print("Running...")
    t_vec_new <- rep(NA,K)
    for(k in 1:K){
      adj_p <- PC_p_list[[k]][selections[[k]]]/weights_list[[k]]
      FDPhatkvec <- FDP_hat_k_vec(PC_p_list = PC_p_list, t_vec = t_vec, k = k,
                                  tk_candidates = sort(adj_p),
                                  pi0_hat_vec = pi0_hat_vec, selections = selections,
                                  weights_list = weights_list, lambdas = lambdas)
      if(sum(FDPhatkvec <= w_vec[k]*q) == 0){
        t_vec_new[k] <- 0
      } else {
        tk_index <- max(which(FDPhatkvec <= w_vec[k]*q))
        t_vec_new[k] <- sort(adj_p)[tk_index]
      }
    }

    t_vec_prev <- t_vec
    t_vec <- t_vec_new
    print(t_vec)
  }

  # Compute the final rejected indices

  R_set <- create_R(PC_p_list = PC_p_list, t_vec = t_vec,
                    selections = selections, weights_list = weights_list,
                    lambdas = lambdas)

  # Output: Vector of rejected indices
  return(R_set)
}


