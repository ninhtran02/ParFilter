#library("partitions", lib.loc = "~/R/lib")
library(partitions)
# Group partitions
#library(partitions)

# COMPUTE ALL NON-EMPTY PARTITIONS OF SIZE k

#' Generate Partitions
#'
#' @param n number of studies
#' @param k maximum number of non-empty sets in the partition
#'
#' @return A list of partitions of \eqn{1,...,n} of up to k non-empty subsets.
#' @export
#'
#' @examples StirlingPartitions(4,3)
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
allStirlingPartions <- function(n,u){
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





