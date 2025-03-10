# COMBINING FUNCTIONS (for testing null PC hypothesis)
# Fisher combining function
fisher_fun <- function(p_vec, u){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level

=======
  
  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  n <- length(p_vec)
  p_subset <- sort(p_vec, decreasing = TRUE)[1:(n - u + 1)]
  len_subset <- length(p_subset)
  PC_p <- 1 - pchisq(q = -2*sum(log(p_subset)), df = 2*len_subset)
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: One PC p-value
  return(PC_p)
}
# Stouffer combining function
stouffer_fun <- function(p_vec, u){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level

=======
  
  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  n <- length(p_vec)
  p_subset <- sort(p_vec, decreasing = TRUE)[1:(n - u + 1)]
  len_subset <- length(p_subset)
  PC_p <- 1 - pnorm(q = sum(qnorm(p = 1 - p_subset))/sqrt(len_subset))
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: One PC p-value
  return(PC_p)
}
# Simes combining function
simes_fun <- function(p_vec, u){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level

  n <- length(p_vec)
  p_subset <- sort(p_vec, decreasing = TRUE)[1:(n - u + 1)]
  PC_p <- min(p.adjust(p = p_subset, method = "BH"))

=======
  
  # Inputs
  ## p_vec: Vector of p-values
  ## u: Replicability level
  
  n <- length(p_vec)
  p_subset <- sort(p_vec, decreasing = TRUE)[1:(n - u + 1)]
  PC_p <- min(p.adjust(p = p_subset, method = "BH"))
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: One PC p-value
  return(PC_p)
}

# DEFAULT SELECTION RULES
default_selections <- function(p_mat, groups, u, method, q, lambdas){
  # Last edited: 22nd October 2024
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Inputs
  ## p_mat: mxn matrix of p-values
  ## groups: List containing the study groupings
  ## u: Replicability level
  ## method: Method for combining p-values: "Fisher", "Stouffer", or "Simes"
  ## q: FDR target
<<<<<<< HEAD
  ## weights: Logical for whether to weight the PC p-values

=======
  ## weights: Logical for whether to weight the PC p-values 
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  K <- length(groups) # Number of groups
  m <- nrow(p_mat) # Number of features
  n <- ncol(p_mat) # Number of studies
  G_sizes <- unlist(lapply(X = groups, FUN = length)) # Group sizes
  min_Gk <- min(G_sizes) # Smallest group size
  u_bar <- max(u - n + min_Gk,1) # Local replicability levels (for selections)
  w_vec <- G_sizes/n # Local error weights
<<<<<<< HEAD

  # Weights for the features
  weights_list <- replicate(K, rep(1,m), simplify = FALSE)

=======
  
  # Weights for the features
  weights_list <- replicate(K, rep(1,m), simplify = FALSE)
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Create a list of local (u_bar/|Gk|) PC p-values
  PC_p_list <- list()
  for(k in 1:K){
    if(u == n){
      u_bar <- G_sizes[k]
      #print(u_bar)
    }
    if(method == "Fisher"){
      PC_p_list[[k]] <- apply(X = cbind(as.matrix(p_mat[,groups[[k]]])), MARGIN = 1, FUN = fisher_fun, u = u_bar)
    }
    if(method == "Stouffer"){
      PC_p_list[[k]] <- apply(X = cbind(as.matrix(p_mat[,groups[[k]]])), MARGIN = 1, FUN = stouffer_fun, u = u_bar)
    }
    if(method == "Simes"){
      PC_p_list[[k]] <- apply(X = cbind(as.matrix(p_mat[,groups[[k]]])), MARGIN = 1, FUN = simes_fun, u = u_bar)
    }
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Threshold the the local PC p-values
  which_thresholded <- list()
  for(k in 1:K){
    which_thresholded[[k]] <- which(PC_p_list[[k]] <= pmin(weights_list[[k]]*w_vec[k]*q,lambdas[k]))
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  selections <- list()
  for(k in 1:K){
    # Intersect thresholded features outside of group k
    selections[[k]] <- Reduce(f = intersect, x = which_thresholded[-k])
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: List of default selections for each group
  return(selections)
}


# COMPUTE REFINED SELECTIONS
refined_selections <- function(p_mat, groups, u_groups_mat, method, q, lambdas){
  # Last edited: 22nd October 2024
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Inputs
  ## p_mat: mxn matrix of p-values
  ## groups: List containing the study groupings
  ## u: Replicability level
  ## method: Method for combining p-values: "Fisher", "Stouffer", or "Simes"
  ## q: FDR target
<<<<<<< HEAD
  ## weights: Logical for whether to weight the PC p-values

=======
  ## weights: Logical for whether to weight the PC p-values 
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  K <- length(groups) # Number of groups
  m <- nrow(p_mat) # Number of features
  n <- ncol(p_mat) # Number of studies
  G_sizes <- unlist(lapply(X = groups, FUN = length)) # Group sizes
  min_Gk <- min(G_sizes) # Smallest group size
  w_vec <- G_sizes/n # Local error weights
<<<<<<< HEAD

  # Weights for the features
  weights_list <- replicate(K, rep(1,m), simplify = FALSE)

=======
  
  # Weights for the features
  weights_list <- replicate(K, rep(1,m), simplify = FALSE)
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Create a list of local (u_bar/|Gk|) PC p-values
  PC_p_list <- list()
  PC_p_list <- list()
  for(k in 1:K){
    if(method == "Fisher"){
      PC_p_vec <- rep(NA, m)
      for(i in 1:m){
        PC_p_vec[i] <- fisher_fun(p_vec = p_mat[i,groups[[k]]], u = u_groups_mat[i,k])
      }
      PC_p_list[[k]] <- PC_p_vec
    }
    if(method == "Stouffer"){
      PC_p_vec <- rep(NA, m)
      for(i in 1:m){
        PC_p_vec[i] <- stouffer_fun(p_vec = p_mat[i,groups[[k]]], u = u_groups_mat[i,k])
      }
      PC_p_list[[k]] <- PC_p_vec
    }
    if(method == "Simes"){
      PC_p_vec <- rep(NA, m)
      for(i in 1:m){
        PC_p_vec[i] <- simes_fun(p_vec = p_mat[i,groups[[k]]], u = u_groups_mat[i,k])
      }
      PC_p_list[[k]] <- PC_p_vec
    }
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Threshold the the local PC p-values
  which_thresholded <- list()
  for(k in 1:K){
    which_thresholded[[k]] <- which(PC_p_list[[k]] <= pmin(weights_list[[k]]*w_vec[k]*q,lambdas[k]))
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  selections <- list()
  for(k in 1:K){
    # Intersect thresholded features outside of group k
    selections[[k]] <- Reduce(f = intersect, x = which_thresholded[-k])
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: List of default selections for each group
  return(selections)
}

# COMPUTE LOCAL REPLICABILITY LEVELS
u_groups_mat_maker <- function(p_mat, X_list, groups, u, lfdr_mat){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Inputs
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## u: Replicability level
<<<<<<< HEAD

  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(groups) # Number of groups

=======
  
  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(groups) # Number of groups
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  if(u == n){
    u_s <- rep(unlist(lapply(groups, length)),m)
    u_groups_mat <- matrix(data = u_s, nrow = m, ncol = K, byrow = TRUE)
    return(u_groups_mat)
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  if(u == K){
    u_groups_mat <- matrix(data = 1, nrow = m, ncol = K)
    return(u_groups_mat)
  }
<<<<<<< HEAD

  # Transform the weights using -log()
  study_weights_mat <- 1 - lfdr_mat

  # Compute the local replicability levels for each feature
  u_groups_mat <- t(apply(X = study_weights_mat, MARGIN = 1,
                          FUN = u_groups_maker, groups = groups, u = u))

=======
  
  # Transform the weights using -log()
  study_weights_mat <- 1 - lfdr_mat
  
  # Compute the local replicability levels for each feature
  u_groups_mat <- t(apply(X = study_weights_mat, MARGIN = 1, 
                          FUN = u_groups_maker, groups = groups, u = u))
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: Matrix of local replicability levels
  return(u_groups_mat)
}

# COMPUTE LOCAL REPLICABILITY LEVELS BASED ON STUDY WEIGHTS
u_groups_maker <- function(study_weights, groups, u){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  ## Inputs
  ## study_weights: Vector of study weights
  ## groups: List of groups
  ## u: Replicability level
<<<<<<< HEAD

  n <- length(study_weights) # Number of studies
  K <- length(groups) # Number of groups

  study_weight_ranks <- rank(-study_weights) # Rank of the weights

=======
  
  n <- length(study_weights) # Number of studies
  K <- length(groups) # Number of groups 
  
  study_weight_ranks <- rank(-study_weights) # Rank of the weights
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  grouped_study_weight_ranks <- list()
  u_groups <- rep(1,K)
  for(k in 1:K){
    group_k <- study_weight_ranks[groups[[k]]]
    group_k <- c(Inf,group_k[-which.min(group_k)])
    grouped_study_weight_ranks[[k]] <- group_k
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  while(sum(u_groups) < u){
    k <- which.min(unlist(lapply(grouped_study_weight_ranks, min)))
    u_groups[k] <- u_groups[k] + 1
    group_k <- grouped_study_weight_ranks[[k]]
    grouped_study_weight_ranks[[k]] <- group_k[-which.min(group_k)]
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: Vector of local replicability levels
  return(u_groups)
}

<<<<<<< HEAD
# COMPUTE APPROXIMATE LFDR
compute_approx_lfdr <- function(p_mat, X_list, groups, selections){
  # Last edited: 23rd October 2024

=======
# COMPUTE APPROXIMATE LFDR 
compute_approx_lfdr <- function(p_mat, X_list, groups, selections){
  # Last edited: 23rd October 2024
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  ## Input
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## selections: List of default selections for each group
<<<<<<< HEAD

  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(selections) # Number of groups

=======
  
  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(selections) # Number of groups
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Make a list of covariates (matrix) for each study
  for(j in 1:n){
    X_list[[j]] <- cbind(X_list[[j]])
  }
<<<<<<< HEAD

  # Create a matrix of estimated values for Pr(H_{ij} = 0 | x_{ij})
  # Get training data from outside the groups
=======
  
  # Create a matrix of estimated values for Pr(H_{ij} = 0 | x_{ij})
  # Get training data from outside the groups 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  training_p_values_list <- list()
  training_covariates_list <- list()
  for(j in 1:n ){
    k_index <- which( unlist(lapply(X = groups, FUN = function(x){ j %in% x })) )
    training_p_values_list[[j]] <- p_mat[-selections[[k_index]],j]
    training_covariates_list[[j]] <- X_list[[j]][-selections[[k_index]],]
  }
<<<<<<< HEAD

  # Train the parameters
  pi0.coef_list <- list()
  k.coef_list <- list()
  for(j in 1:n){
    camt_obj <- CAMT::camt.fdr(pvals = training_p_values_list[[j]],
=======
  
  # Train the parameters 
  pi0.coef_list <- list()
  k.coef_list <- list()
  for(j in 1:n){
    camt_obj <- CAMT::camt.fdr(pvals = training_p_values_list[[j]], 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
                               pi0.var = training_covariates_list[[j]],
                               f1.var = training_covariates_list[[j]])
    pi0.coef_list[[j]] <- 1.5*camt_obj$pi0.coef
    k.coef_list[[j]]   <- 1.0*camt_obj$k.coef
  }
<<<<<<< HEAD


=======
  
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
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
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  p_mat_for_weights <- -(k_est_mat*pi0_est_mat - 2*k_est_mat + 2)/(2*k_est_mat - 4)
  hp_mat <- (1 - k_est_mat)*p_mat_for_weights^(-k_est_mat)
  hp_mat <- (1 - k_est_mat)*pmin(p_mat_for_weights^(-k_est_mat),100)
  lfdr_mat <- pi0_est_mat/( (1 - pi0_est_mat)*hp_mat + pi0_est_mat )
<<<<<<< HEAD

  #lfdr_mat <- 1 - k_est_mat
  #lfdr_mat <- pi0_est_mat

=======
  
  #lfdr_mat <- 1 - k_est_mat
  #lfdr_mat <- pi0_est_mat
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  pi0.coef <- c()
  k.coef <- c()
  for(k in 1:K){
    pi0.coef <- c(pi0.coef, pi0.coef_list[[k]][-1])
    k.coef <-   c(k.coef,   k.coef_list[[k]][-1]  )
  }
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
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
<<<<<<< HEAD

  print("PROTOTYPE K")

=======
  
  print("PROTOTYPE K")
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  ## Input
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## u_groups_mat: Matrix of local replicability levels
  ## selections: List of refined selections for each group
<<<<<<< HEAD
  K <- length(groups)
  weights_list <- list() # List of weights for selected features for each group
=======
  
  weights_list <- list() # List of weights for selected features for each group 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571

  for(k in 1:K){
    est_ulfdr_vec <- rep(0, length(selections[[k]]))
    all_combinations <- expand.grid(rep(list(c(0, 1)), length(groups[[k]]) ))
    for(i in 1:length(selections[[k]])){
      valid_combinations <- cbind(all_combinations[rowSums(all_combinations) < u_groups_mat[ selections[[k]][i] ,k], ])
<<<<<<< HEAD

=======
      
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
      for( l in 1:nrow(valid_combinations) ){
        est_ulfdr_part_zeros <- prod(lfdr_mat[ selections[[k]][i],    groups[[k]][which(valid_combinations[l,] == 0)]  ]    )
        est_ulfdr_part_ones <- prod(1- lfdr_mat[ selections[[k]][i] , groups[[k]][which(valid_combinations[l,] == 1)]  ]    )
        est_ulfdr_vec[i] <- est_ulfdr_vec[i] + est_ulfdr_part_zeros*est_ulfdr_part_ones
      }
    }
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
    # Scale the weights
    if(normalize){
      feature_weights <- 1 - est_ulfdr_vec
      feature_weights <- (feature_weights/sum(feature_weights))*length(selections[[k]])
      weights_list[[k]] <- feature_weights
    } else {
      weights_list[[k]] <- 1 - est_ulfdr_vec
    }
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: List of weights for selected features in each group
  return(weights_list)
}

# COMPUTE WEIGHTS FOR FDR CONTROL
cross_weights_maker <- function(lfdr_mat, groups, u_groups_mat, selections,
                          normalize = TRUE){
  # Last edited: 23rd October 2024
<<<<<<< HEAD

  print("PROTOTYPE K")

=======
  
  print("PROTOTYPE K")
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  ## Input
  ## p_mat: mxn matrix of p-values
  ## X_list: List of covariates for each study (length(X_list) should be n)
  ## groups: List containing the study groupings
  ## u_groups_mat: Matrix of local replicability levels
  ## selections: List of refined selections for each group
<<<<<<< HEAD

  n <- ncol(cbind(lfdr_mat))

  weights_list <- list() # List of weights for selected features for each group

=======
  
  n <- ncol(cbind(lfdr_mat))
  
  weights_list <- list() # List of weights for selected features for each group 
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  for(k in 1:K){
    est_ulfdr_vec <- rep(0, length(selections[[k]]))
    all_combinations <- expand.grid(rep(list(c(0, 1)), n - length(groups[[k]]) ))
    oosg <- unlist(groups[-k])
    for(i in 1:length(selections[[k]])){
      valid_combinations <- cbind(all_combinations[rowSums(all_combinations) < u - u_groups_mat[ selections[[k]][i] ,k], ])
<<<<<<< HEAD

=======
      
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
      for( l in 1:nrow(valid_combinations) ){
        est_ulfdr_part_zeros <- prod(lfdr_mat[ selections[[k]][i],    oosg[which(valid_combinations[l,] == 0)]  ]    )
        est_ulfdr_part_ones <- prod(1- lfdr_mat[ selections[[k]][i] , oosg[which(valid_combinations[l,] == 1)]  ]    )
        est_ulfdr_vec[i] <- est_ulfdr_vec[i] + est_ulfdr_part_zeros*est_ulfdr_part_ones
      }
    }
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
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
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: List of weights for selected features in each group
  return(weights_list)
}

# COMPUTE PC P-VALUES FOR EACH GROUP
PC_p_maker <- function(p_mat, groups, u_groups_mat, method){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  ## Inputs
  ## p_mat: mxn matrix of p-values
  ## groups: List containing the study groupings
  ## u_groups_mat: Matrix of local replicability levels
  ## method: Method for combining p-values: "Fisher", "Stouffer", or "Simes"
<<<<<<< HEAD

  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(groups) # Number of groups

=======
  
  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  K <- length(groups) # Number of groups
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
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
<<<<<<< HEAD

  # Output: List of PC p-values for each group
=======
  
  # Output: List of PC p-values for each group 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  return(PC_p_list)
}

# COMPUTE (WEIGHTED) STOREY NULL PROPORTIONE ESTIMATE (2004)
pi0_hat_maker <- function(PC_p_list,
                    selections, weights_list, lambdas){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  ## Inputs:
  ## PC_p_list: List of PC p-values for each group
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters

=======
  
  ## Inputs:
  ## PC_p_list: List of PC p-values for each group 
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  K <- length(selections) # Number of groups
  pi0_hat_vec <- rep(NA,K)
  for(k in 1:K){
    num_sum_index <- which(PC_p_list[[k]][selections[[k]]] > lambdas[k])
    numerator <- max(weights_list[[k]]) + sum(weights_list[[k]][num_sum_index])
    denominator <- length(selections[[k]])*(1 - lambdas[k])
    pi0_hat_vec[k] <- numerator/denominator
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: Vector of (post-selection) pi0 estimates
  return(pi0_hat_vec)
}

# COMPUTE (WEIGHTED) STOREY NULL PROPORTIONE ESTIMATE (2002)
pi0_hat_maker_storey_2002 <- function(PC_p_list,
                          selections, weights_list, lambdas){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  ## Inputs:
  ## PC_p_list: List of PC p-values for each group
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters

=======
  
  ## Inputs:
  ## PC_p_list: List of PC p-values for each group 
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  K <- length(selections) # Number of groups
  pi0_hat_vec <- rep(NA,K)
  for(k in 1:K){
    num_sum_index <- which(PC_p_list[[k]][selections[[k]]] > lambdas[k])
    numerator <- sum(weights_list[[k]][num_sum_index])
    denominator <- length(selections[[k]])*(1 - lambdas[k])
    pi0_hat_vec[k] <- numerator/denominator
  }
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: Vector of (post-selection) pi0 estimates
  return(pi0_hat_vec)
}

# BOOSTRAP BEST VALUES FOR LAMBDAS
auto_lambda <-  function(PC_p_list,
                         selections, weights_list, q, w_vec, B = 1000){
  K <- length(selections)
  optimal_lambdas <- rep(NA, K)
  m <- length(PC_p_list[[1]])
  for(k in 1:K){
    lambda_options <- seq(0.05, 0.95, 0.05)
    pi0_hat_storey_2002_vec <- rep(NA, length(lambda_options))
    for(i in 1:length(lambda_options)){
<<<<<<< HEAD
      #pi0_hat_storey_2002_vec[i] <- pi0_hat_maker_storey_2002(PC_p_list = PC_p_list[k],
=======
      #pi0_hat_storey_2002_vec[i] <- pi0_hat_maker_storey_2002(PC_p_list = PC_p_list[k], 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
      #                                                        selections = selections[k],
      #                                                        weights_list = weights_list[k],
      #                                                        lambdas = lambda_options[i])
      unselected_indices <- setdiff(1:m,selections[[k]])
      pi0_hat_storey_2002_vec[i] <- sum(PC_p_list[[k]][unselected_indices] > lambda_options[i])/((1-lambda_options[i])*length(unselected_indices))
    }
    min_pi0_hat_storey_2002 <- min(pi0_hat_storey_2002_vec)
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
    MSE_vec <- rep(NA, length(lambda_options))
    for(i in 1:length(lambda_options)){
      MSE <- 0
      for(b in 1:B){
        boostrap_p <- sample(x = PC_p_list[[k]][unselected_indices], size = length(unselected_indices), replace = TRUE)
        pi0_hat <- (1 + sum(boostrap_p > lambda_options[i]))/((1-lambda_options[i])*length(unselected_indices))
<<<<<<< HEAD

        #pi0_hat <- pi0_hat_maker(PC_p_list = PC_p_list[k],
=======
        
        #pi0_hat <- pi0_hat_maker(PC_p_list = PC_p_list[k], 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
        #                         selections = selections[k],
        #                         weights_list = weights_list[k],
        #                         lambdas = lambda_options[i])
        MSE <- MSE + (1/B)*(pi0_hat - min_pi0_hat_storey_2002)^2
      }
      MSE_vec[i] <- MSE
    }
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
    optimal_lambdas[k] <- lambda_options[which.min(MSE_vec)]
  }
  return(optimal_lambdas)
}


# COMPUTE ParFilter REJECTION SET
create_R <- function(PC_p_list, t_vec, selections,
                     weights_list, lambdas){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  ## Inputs
  ## PC_p_list: List of PC p-values for each group
=======
  
  ## Inputs
  ## PC_p_list: List of PC p-values for each group 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  ## t_vec: Vector of rejection thresholds for each group
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
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
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: Vector of rejected indices
  return(R_set)
}

# Compute a vector of FDP_hat_k values based on a
# vector of candidate values for the group k threshold
FDP_hat_k_vec <- function(PC_p_list, t_vec, k, tk_candidates,
                          pi0_hat_vec, selections,
                          weights_list, lambdas){
<<<<<<< HEAD

  # Last edited: 3rd October 2024

  ## Inputs
  ## PC_p_list: List of PC p-values for each group
=======
  
  # Last edited: 3rd October 2024
  
  ## Inputs
  ## PC_p_list: List of PC p-values for each group 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  ## t_vec: Vector of rejection thresholds for each group
  ## k: Group of interest
  ## tk_candidates: Vector of values for tk
  ## pi0_hat_vec: Vector of (post-selection) pi0 estimates
  ## selections: List of refined selected features for each group
  ## weights_list: List of weights for selected features in each group
  ## lambdas: Vector of hyperparameters
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  R_sizes <- rep(NA, length(tk_candidates)) # Vector to store rejection sizes
  for(i in 1:length(tk_candidates) ){
    t_vec_temp <- t_vec
    t_vec_temp[k] <- tk_candidates[i] # Replace current tk with candidate
    R_set <- create_R(PC_p_list = PC_p_list, t_vec = t_vec_temp,
                      selections = selections, weights_list = weights_list,
                      lambdas = lambdas)
    R_sizes[i] <- length(R_set)
  }
<<<<<<< HEAD

  FDP_hat_vec <- length(selections[[k]])*pi0_hat_vec[k]*tk_candidates/pmax(R_sizes,1)

=======
  
  FDP_hat_vec <- length(selections[[k]])*pi0_hat_vec[k]*tk_candidates/pmax(R_sizes,1)
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: Vector of FDP_hat values for group k given candidate tk thresholds
  return(FDP_hat_vec)
}

# COMPUTE GROUPS
groups_maker <- function(u,n,K){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  group_index <- cut(1:n, breaks = K, labels = FALSE)
  groups <- split(1:n, group_index)

=======
  
  group_index <- cut(1:n, breaks = K, labels = FALSE)
  groups <- split(1:n, group_index)
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: List of study groupings
  return(groups)
}

# COMPUTE LOCAL ERROR WEIGHTS
w_vec_maker <- function(groups){
  # Last edited: 3rd October 2024
<<<<<<< HEAD

  w_vec <- unlist(lapply(X = groups, FUN = length))/length(unlist(groups))

=======
  
  w_vec <- unlist(lapply(X = groups, FUN = length))/length(unlist(groups))
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Output: Vector of local error weights
  return(w_vec)
}

# ParFilter ALGORITHM
<<<<<<< HEAD
#' ParFilter
#'
#' @param p_mat mxn matrix of p-values.
#' @param X_list list of length n, containing the covariates for each study.
#' @param u replicability threshold.
#' @param q FDR target.
#' @param K number of groups. ParFilter_FDR will automatically partition the n studies in two K groups of approximately equal sizes.
#' @param method combining function for creating the local GBHPC p-values. Can be either: "Fisher", "Stouffer", or "Simes".
#' @param adaptive logical indicating whether to use adaptive null proportion estimates or not.
#' @param lambdas set as TRUE if the p-values are dependent within studies, otherwise leave it as FALSE.
#' @param cross_weights numeric of tuning parameters for the null proportion estimates.

#'
#' @return returns a numeric of features indices to be rejected.
#' @export
#'
#' @examples
=======
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
ParFilter_FDR <- function(p_mat, X_list, u, q, K, method,
                          adaptive = TRUE, lambdas = rep(0.5,K),
                          cross_weights = FALSE){
  # Last edited: 22nd October 2024
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
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
<<<<<<< HEAD

  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features

  # Compute the study groupings and local error_weights
  groups <- groups_maker(u = u, n = n, K = K)
  w_vec <- w_vec_maker(groups)

=======
  
  n <- ncol(p_mat) # Number of studies
  m <- nrow(p_mat) # Number of features
  
  # Compute the study groupings and local error_weights
  groups <- groups_maker(u = u, n = n, K = K)
  w_vec <- w_vec_maker(groups)
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Form the selection rules
  if(u == n){
    selections <- default_selections(p_mat = p_mat, groups = groups,
                                     u = u, method = method, q = q,
                                     lambdas = lambdas)
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  } else {
    selections <- default_selections(p_mat = p_mat, groups = groups,
                                     u = u, method = method, q = q,
                                     lambdas = lambdas)
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
    #lfdr_obj <- compute_approx_lfdr(p_mat = p_mat, X_list = X_list,
    #                                groups = groups,
    #                                selections = selections)
    #lfdr_mat <- lfdr_obj$lfdr_mat
    #u_groups_mat <- u_groups_mat_maker(p_mat = p_mat, X_list = X_list,
    #                                   groups = groups, u = u,
    #                                   lfdr_mat = lfdr_mat)
    #selections <- refined_selections(p_mat = p_mat, groups = groups,
    #                                 u_groups_mat = u_groups_mat, method = method, q = q,
    #                                 lambdas = lambdas)
  }
<<<<<<< HEAD


  if(0 %in% lapply(selections, length)){
    return(c())
  }

=======
  
  
  if(0 %in% lapply(selections, length)){
    return(c())
  }
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  lfdr_obj <- compute_approx_lfdr(p_mat = p_mat, X_list = X_list,
                                  groups = groups,
                                  selections = selections)
  lfdr_mat <- lfdr_obj$lfdr_mat
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Make u_groups_mat
  u_groups_mat <- u_groups_mat_maker(p_mat = p_mat, X_list = X_list,
                                     groups = groups, u = u,
                                     lfdr_mat = lfdr_mat)
<<<<<<< HEAD

=======
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Form the weights
  if(cross_weights == "naive"){
    weights_list <- list()
    for(k in 1:K){
      weights_list[[k]] <- rep(1,length(selections[[k]]))
    }
  }
  if(cross_weights == FALSE){
<<<<<<< HEAD
    weights_list <- weights_maker(lfdr_mat = lfdr_mat, groups = groups,
=======
    weights_list <- weights_maker(lfdr_mat = lfdr_mat, groups = groups, 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
                                  u_groups_mat = u_groups_mat,
                                  selections = selections)
  }
  if(cross_weights == TRUE){
<<<<<<< HEAD
    weights_list <- cross_weights_maker(lfdr_mat = lfdr_mat, groups = groups,
                                  u_groups_mat = u_groups_mat,
                                  selections = selections)
  }


  # Form the PC p-values
  PC_p_list <- PC_p_maker(p_mat = p_mat, groups = groups,
                          u_groups_mat = u_groups_mat, method = method)

  #print(lambdas)

=======
    weights_list <- cross_weights_maker(lfdr_mat = lfdr_mat, groups = groups, 
                                  u_groups_mat = u_groups_mat,
                                  selections = selections)
  }
  
  
  # Form the PC p-values
  PC_p_list <- PC_p_maker(p_mat = p_mat, groups = groups, 
                          u_groups_mat = u_groups_mat, method = method)

  #print(lambdas)
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  # Form the pi0_hats
  if(adaptive){
    #for( k in 1:K ){
    #  lambdas[k] <- 1 - min(mean(lfdr_obj$pi0_est_mat[selections[[k]],k]*weights_list[[k]]),0.95)
    #}
    #print(lambdas)
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
    pi0_hat_vec <- pi0_hat_maker(PC_p_list = PC_p_list, selections = selections,
                                 weights_list = weights_list, lambdas = lambdas)
  } else {
    pi0_hat_vec <- rep(1,K)
  }
<<<<<<< HEAD

  # Find the rejection thresholds
  t_vec <- rep( m, K )
  t_vec_prev <- rep(0,K)

=======
  
  # Find the rejection thresholds
  t_vec <- rep( m, K )
  t_vec_prev <- rep(0,K)
  
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
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
<<<<<<< HEAD

=======
    
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
    t_vec_prev <- t_vec
    t_vec <- t_vec_new
    print(t_vec)
  }
<<<<<<< HEAD

  # Compute the final rejected indices

  R_set <- create_R(PC_p_list = PC_p_list, t_vec = t_vec,
                    selections = selections, weights_list = weights_list,
                    lambdas = lambdas)

  # Output: Vector of rejected indices
=======
  
  # Compute the final rejected indices
  
  R_set <- create_R(PC_p_list = PC_p_list, t_vec = t_vec,
                    selections = selections, weights_list = weights_list,
                    lambdas = lambdas)
  
  # Output: Vector of rejected indices 
>>>>>>> b79beae35f69c05bb92e0182f7fddc6e96c94571
  return(R_set)
}


