# Set the library path

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


library("partitions", lib.loc = "~/R/lib")
library("adaFilter")

source("auto.R")
source("RepFilterIndependence.R")
source("stirlingfunctions.R")
source("simulation_functions.R")

cmd_args <- commandArgs(TRUE)
rho_index = as.numeric(cmd_args[1])
mu_index = as.numeric(cmd_args[2])
setup_index = as.numeric(cmd_args[3])

rho_options <- c(-0.8,0,0.8)
u_n_settings <- list(c(2,2),c(3,4),c(4,4),
                     c(4,6),c(5,6),c(6,6),
                     c(6,8),c(7,8),c(8,8))
#mu_options <- c(2,2.5,3,3.5,4)
mu_options <- c(2,2.3,2.6,2.9,3.2,3.5)

rho <- rho_options[rho_index]
setup <- u_n_settings[[setup_index]]
mu <- mu_options[mu_index]

methods <- c("PF", "MaxPF", "APF", "MaxAPF", "AutoAPF", "BH-Simes",
             "BH-Stouffer", "BH-Fisher", "AdaFilter")

if(rho != 0){
  methods <- c(methods, "BY-Simes", "BY-Stouffer", "BY-Fisher", "ad-PF", "prds-PF" )
}




### SETTINGS
# Number of rows
m <- 5000
# Number of simulations
nsim <- 200
alpha <- 0.05

set.seed(832649)

u <- setup[1]
n <- setup[2]
group_options <- allStirlingPartitions(n = n, u = u)[-1]
#group_options <- SafePartitions(n = n, u = u)
#if(n == 8){
#  group_options <- group_options[which(unlist(lapply(group_options, length)) == 2)]
#}

FDRrep <- rep(0,length(methods))
names(FDRrep) <- methods

TPRrep <- rep(0,length(methods))
names(TPRrep) <- methods

FDRstudy <- matrix(0,nrow = n, ncol = length(methods))
colnames(FDRstudy) <- methods

TPRstudy <- matrix(0,nrow = n, ncol = length(methods))
colnames(TPRstudy) <- methods

print(paste("rho",rho,"mu",mu,"setup",u,n,sep = ""))

for(iter in 1:nsim){
  print(iter)
  rep_rejection_list <- vector("list", length(methods))
  study_rejection_list <- vector("list", length(methods))
  names(rep_rejection_list) <- methods
  data_list <- generate_pvalues(m = m, n = n, u = u, mu = mu, rho = rho)
  p <- data_list$p
  class(p) <- "matrix"

  for(method in methods){
    if(method == "PF"){
      obj <- parfilter(p = p, error_targets = c(rep(alpha/5, n),alpha),
                      u = u, groups = NULL, group_options = group_options,
                      selections = NULL, u_groups = NULL,
                      adaptive = FALSE, lambda = NULL,
                      auto = FALSE, omega = 0.5, w = NULL,
                      partition = "minimum", method = "Stouffer")
      rep_rejection_list[[method]] <- obj$Replicability_Rejections
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "MaxPF"){
      obj <- parfilter(p = p, error_targets = c(rep(alpha/5, n),alpha),
                       u = u, groups = NULL, group_options = group_options,
                       selections = NULL, u_groups = NULL,
                       adaptive = FALSE, lambda = NULL,
                       auto = FALSE, omega = 0.5, w = NULL,
                       partition = "maximum", method = "Stouffer")
      rep_rejection_list[[method]] <- obj$Replicability_Rejections
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "APF"){
      obj <- parfilter(p = p, error_targets = c(rep(alpha/5, n),alpha),
                       u = u, groups = NULL, group_options = group_options,
                       selections = NULL, u_groups = NULL,
                       adaptive = TRUE, lambda = NULL,
                       auto = FALSE, omega = 0.5, w = NULL,
                       partition = "minimum", method = "Stouffer")
      rep_rejection_list[[method]] <- obj$Replicability_Rejections
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "MaxAPF"){
      obj <- parfilter(p = p, error_targets = c(rep(alpha/5, n),alpha),
                       u = u, groups = NULL, group_options = group_options,
                       selections = NULL, u_groups = NULL,
                       adaptive = TRUE, lambda = NULL,
                       auto = FALSE, omega = 0.5, w = NULL,
                       partition = "maximum", method = "Stouffer")
      rep_rejection_list[[method]] <- obj$Replicability_Rejections
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "AutoAPF"){
      #if(length(group_options) > 100){
      #  sample_groups <- sample(1:length(group_options),
      #                          size = 100, replace = FALSE)
      #} else {
      #  sample_groups <- 1:length(group_options)
      #}
      groups <- group_options[sample(x = length(group_options), size = 1)][[1]]
      #obj <- parfilter(p = p, error_targets = rep(alpha, n + 1),
      #                 u = u, groups = NULL, group_options = group_options[sample_groups],
      #                 selections = NULL, u_groups = NULL,
      #                 adaptive = TRUE, lambda = NULL,
      #                 auto = TRUE, omega = 0.75, w = NULL,
      #                 partition = NULL, method = "Stouffer")
      obj <- parfilter(p = p, error_targets = c(rep(alpha/5, n),alpha),
                       u = u, groups = groups, group_options = group_options,
                       selections = NULL, u_groups = NULL,
                       adaptive = TRUE, lambda = NULL,
                       auto = FALSE, omega = 0.5, w = NULL,
                       partition = NULL, method = "Stouffer")
      rep_rejection_list[[method]] <- obj$Replicability_Rejections
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "BH-Simes"){
      p_list <- split(c(t(p)), rep(1:nrow(p), each = ncol(p)))
      u_list <- split(rep(u,m), 1:m)
      method_list <- split(rep("Simes",m), 1:m)
      simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
      simes_pu <- unname(simes_pu)
      rep_rejection_list[[method]] <- which(p.adjust(p = simes_pu,
                                                     method = "BH") <= alpha)
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "BH-Stouffer"){
      p_list <- split(c(t(p)), rep(1:nrow(p), each = ncol(p)))
      u_list <- split(rep(u,m), 1:m)
      method_list <- split(rep("Stouffer",m), 1:m)
      simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
      simes_pu <- unname(simes_pu)
      rep_rejection_list[[method]] <- which(p.adjust(p = simes_pu,
                                                     method = "BH") <= alpha)
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "BH-Fisher"){
      p_list <- split(c(t(p)), rep(1:nrow(p), each = ncol(p)))
      u_list <- split(rep(u,m), 1:m)
      method_list <- split(rep("Fisher",m), 1:m)
      simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
      simes_pu <- unname(simes_pu)
      rep_rejection_list[[method]] <- which(p.adjust(p = simes_pu,
                                                     method = "BH") <= alpha)
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "AdaFilter"){
      obj <- adaFilter(p.matrix = p, r = u, type.I.err = "FDR",
                       alpha = alpha/sum(1/1:m))
      rep_rejection_list[[method]] <- which(obj$decision == 1)
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "BY-Simes"){
      p_list <- split(c(t(p)), rep(1:nrow(p), each = ncol(p)))
      u_list <- split(rep(u,m), 1:m)
      method_list <- split(rep("Simes",m), 1:m)
      simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
      simes_pu <- unname(simes_pu)
      rep_rejection_list[[method]] <- which(p.adjust(p = simes_pu,
                                                     method = "BY") <= alpha)
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "BY-Stouffer"){
      p_list <- split(c(t(p)), rep(1:nrow(p), each = ncol(p)))
      u_list <- split(rep(u,m), 1:m)
      method_list <- split(rep("Stouffer",m), 1:m)
      simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
      simes_pu <- unname(simes_pu)
      rep_rejection_list[[method]] <- which(p.adjust(p = simes_pu,
                                                     method = "BY") <= alpha)
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "BY-Fisher"){
      p_list <- split(c(t(p)), rep(1:nrow(p), each = ncol(p)))
      u_list <- split(rep(u,m), 1:m)
      method_list <- split(rep("Fisher",m), 1:m)
      simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
      simes_pu <- unname(simes_pu)
      rep_rejection_list[[method]] <- which(p.adjust(p = simes_pu,
                                                     method = "BY") <= alpha)
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "ad-PF"){
      obj <- parfilter(p = p, error_targets = rep(alpha/sum(1/1:m), n + 1),
                       u = u, groups = NULL, group_options = group_options,
                       selections = NULL, u_groups = NULL,
                       adaptive = FALSE, lambda = NULL,
                       auto = FALSE, omega = 0.5, w = NULL,
                       partition = "minimum", method = "Stouffer")
      rep_rejection_list[[method]] <- obj$Replicability_Rejections
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

    if(method == "prds-PF"){
      obj <- parfilter(p = p, error_targets = c(rep(alpha/5, n),alpha),
                       u = u, groups = NULL, group_options = group_options,
                       selections = "PRDS", u_groups = NULL,
                       adaptive = FALSE, lambda = NULL,
                       auto = FALSE, omega = 0.5, w = NULL,
                       partition = "minimum", method = "Stouffer")
      rep_rejection_list[[method]] <- obj$Replicability_Rejections
      study_rejection_list[[method]] <- obj$Studywise_Rejections
    }

  }

  fdprep_numeric <- unlist(lapply(X = rep_rejection_list, FUN = compute_FDP,
         true_nulls = setdiff(1:m,data_list$rep_nonnull_indices)))

  tpprep_numeric <- unlist(lapply(X = rep_rejection_list, FUN = compute_TPP,
         false_nulls = data_list$rep_nonnull_indices))

  fdrstudy <- matrix(0,nrow = n, ncol = length(methods))
  colnames(fdrstudy) <- methods
  tprstudy <- matrix(0,nrow = n, ncol = length(methods))
  colnames(tprstudy) <- methods

  for(method in methods){
    for(j in 1:n){
      fdrstudy[j, method] <- compute_FDP(rejections = study_rejection_list[[method]][[j]],
                                         true_nulls  = setdiff(1:m,data_list$study_nonnull_indices[[j]]))

      tprstudy[j, method] <- compute_TPP(rejections = study_rejection_list[[method]][[j]],
                                         false_nulls = data_list$study_nonnull_indices[[j]])
    }
  }

  FDRrep <- FDRrep + fdprep_numeric/nsim
  TPRrep <- TPRrep + tpprep_numeric/nsim

  FDRstudy <- FDRstudy + fdrstudy/nsim
  TPRstudy <- TPRstudy + tprstudy/nsim
}

if(rho == 0){
  file.name <- paste("SavedData/Independence/","mu",mu,"setup",u,n,".RD",sep = "")
} else if(rho > 0){
  file.name <- paste("SavedData/PositiveDependence/","mu",mu,"setup",u,n,".RD",sep = "")
} else if(rho < 0){
  file.name <- paste("SavedData/NegativeDependence/","mu",mu,"setup",u,n,".RD",sep = "")
}
save(FDRrep, TPRrep, FDRstudy, TPRstudy, file = file.name)

