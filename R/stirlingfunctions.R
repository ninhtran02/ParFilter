library("partitions", lib.loc = "~/R/lib")

# Group partitions
#library(partitions)

# COMPUTE ALL NON-EMPTY PARTITIONS OF SIZE k
StirlingPartitions <- function(n,k){
  # INPUTS
  # n is the number of studies
  # k is the number of non-empty sets in the partition
  
  # Last modified: Ninh Tran 27 Feb 2024
  
  groups_candidates <- parts(n)
  groups_candidates <- cbind(groups_candidates[1:k,which(colSums(groups_candidates == 0) == n - k)])
  groups_candidates <- as.list(data.frame(groups_candidates))
  group_options <- list()
  
  # This a technical step if k = 1
  if(length(groups_candidates) == 1 & length(groups_candidates[[1]]) == 1){ return(list(1:n))}
  
  for(g in groups_candidates){
    sp <- setparts(x = g)
    for(l in 1:dim(sp)[2]){
      group_options <- c(group_options, list(split(x = 1:n, f = sp[,l])))
    }
  }
  
  return(group_options)
}

# COMPUTE ALL NON-EMPTY PARTITIONS
allStirlingPartitions <- function(n,u){
  # INPUTS
  # n is the number of studies
  # u is the maximum number of non-empty sets in the partition
  
  # Last modified: Ninh Tran 27 Feb 2024
  
  group_options <- list()
  for(k in 1:u){
    group_options <- c(group_options, StirlingPartitions(n = n, k = k))
  }
  
  return(group_options)
}

SafePartitions <- function(n,u){
  group_options <- allStirlingPartitions(n = n, u = u)
  safe_indices <- c()
  for(i in 1:length(group_options)){
    len_vec <- sapply(X = group_options[[i]], FUN = length)
    if( cumsum(len_vec[-length(len_vec)])[length(len_vec)-1] < u ){
      safe_indices <- c(safe_indices, i)
    }
  }
  return(group_options[safe_indices])
}




