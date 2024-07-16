#devtools::install_github("ninhtran02/Parfilter")
#devtools::install("jingshuw/adaFilter")
library(ParFilter)
library(adaFilter)

# Load the base p-values
load("DMD.pvalues.Rd")
p.mat <- DMD.pvalues
n <- dim(p.mat)[2]
m <- dim(p.mat)[1]

methods <- c("PF", "MaxPF", "APF", "MaxAPF", "BH-Simes",
             "BH-Stouffer", "BH-Fisher", "AdaFilter")
alpha_options <- c(0.01, 0.05, 0.10, 0.15, 0.20)
repRejection_df <- data.frame(data.frame(matrix(ncol = 4, nrow = 0)))
colnames(repRejection_df) <- c("rejections", "u", "method", "alpha")
studyRejectiondf <- data.frame(data.frame(matrix(ncol = 5, nrow = 0)))
colnames(studyRejectiondf) <- c("rejections", "study", "u", "method", "alpha")

for(u in 2:n){
  print(u)
  for(method in methods){
    for(alpha in alpha_options){
      print(method)
      print(alpha)
      if(method == "PF"){
        groups <- list(c(3,4),c(1,2))
        obj <- parfilter(p = p.mat, error_targets = c(rep(alpha/5, n),alpha),
                         u = u, groups = groups, group_options = group_options,
                         selections = NULL, u_groups = NULL,
                         adaptive = FALSE, lambda = NULL,
                         auto = FALSE, omega = 0.75, w = NULL,
                         partition = NULL, method = "Stouffer")
        data.row <- c(length(obj$Replicability_Rejections), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
        for(j in 1:n){
          data.row <- c(length(obj$Studywise_Rejections[[j]]), j, u, method, alpha)
          studyRejectiondf <- rbind.data.frame(studyRejectiondf, data.row)
        }
      }
      
      if(method == "MaxPF"){
        if(u == 2){
          groups <- list(c(3,4),c(1,2))
        }
        if(u == 3){
          groups <- list(c(1,2),3,4)
        }
        if(u == 4){
          groups <- list(1,2,3,4)
        }
        obj <- parfilter(p = p.mat, error_targets = c(rep(alpha/5, n),alpha),
                         u = u, groups = groups, group_options = group_options,
                         selections = NULL, u_groups = NULL,
                         adaptive = FALSE, lambda = NULL,
                         auto = FALSE, omega = 0.75, w = NULL,
                         partition = NULL, method = "Stouffer")
        data.row <- c(length(obj$Replicability_Rejections), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
        for(j in 1:n){
          data.row <- c(length(obj$Studywise_Rejections[[j]]), j, u, method, alpha)
          studyRejectiondf <- rbind.data.frame(studyRejectiondf, data.row)
        }
      }
      
      if(method == "APF"){
        groups <- list(c(3,4),c(1,2))
        obj <- parfilter(p = p.mat, error_targets = c(rep(alpha/5, n),alpha),
                         u = u, groups = groups, group_options = group_options,
                         selections = NULL, u_groups = NULL,
                         adaptive = TRUE, lambda = NULL,
                         auto = FALSE, omega = 0.75, w = NULL,
                         partition = NULL, method = "Stouffer")
        data.row <- c(length(obj$Replicability_Rejections), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
        for(j in 1:n){
          data.row <- c(length(obj$Studywise_Rejections[[j]]), j, u, method, alpha)
          studyRejectiondf <- rbind.data.frame(studyRejectiondf, data.row)
        }
      }
      
      if(method == "MaxAPF"){
        if(u == 2){
          groups <- list(c(3,4),c(1,2))
        }
        if(u == 3){
          groups <- list(c(1,2),3,4)
        }
        if(u == 4){
          groups <- list(1,2,3,4)
        }
        obj <- parfilter(p = p.mat, error_targets = c(rep(alpha/5, n),alpha),
                         u = u, groups = groups, group_options = group_options,
                         selections = NULL, u_groups = NULL,
                         adaptive = TRUE, lambda = NULL,
                         auto = FALSE, omega = 0.75, w = NULL,
                         partition = NULL, method = "Stouffer")
        data.row <- c(length(obj$Replicability_Rejections), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
        for(j in 1:n){
          data.row <- c(length(obj$Studywise_Rejections[[j]]), j, u, method, alpha)
          studyRejectiondf <- rbind.data.frame(studyRejectiondf, data.row)
        }
      }
      
      if(method == "AutoAPF"){
        obj <- parfilter(p = p.mat, error_targets = c(rep(alpha/5, n),alpha),
                         u = u, groups = NULL, group_options = group_options,
                         selections = NULL, u_groups = NULL,
                         adaptive = TRUE, lambda = NULL,
                         auto = TRUE, omega = 0.75, w = NULL,
                         partition = NULL, method = "Stouffer")
        data.row <- c(length(obj$Replicability_Rejections), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
        for(j in 1:n){
          data.row <- c(length(obj$Studywise_Rejections[[j]]), j, u, method, alpha)
          studyRejectiondf <- rbind.data.frame(studyRejectiondf, data.row)
        }
      }
      
      if(method == "BH-Simes"){
        p_list <- split(c(t(p.mat)), rep(1:nrow(p.mat), each = ncol(p.mat)))
        u_list <- split(rep(u,m), 1:m)
        method_list <- split(rep("Simes",m), 1:m)
        simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
        simes_pu <- unname(simes_pu)
        data.row <- c(length(which(p.adjust(p = simes_pu,
                                            method = "BH") <= alpha)), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
      }
      
      if(method == "BH-Fisher"){
        p_list <- split(c(t(p.mat)), rep(1:nrow(p.mat), each = ncol(p.mat)))
        u_list <- split(rep(u,m), 1:m)
        method_list <- split(rep("Fisher",m), 1:m)
        simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
        simes_pu <- unname(simes_pu)
        data.row <- c(length(which(p.adjust(p = simes_pu,
                                            method = "BH") <= alpha)), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
      }
      
      if(method == "BH-Stouffer"){
        p_list <- split(c(t(p.mat)), rep(1:nrow(p.mat), each = ncol(p.mat)))
        u_list <- split(rep(u,m), 1:m)
        method_list <- split(rep("Stouffer",m), 1:m)
        simes_pu <- mapply(PC_u, p = p_list, u = u_list, method = method_list)
        simes_pu <- unname(simes_pu)
        data.row <- c(length(which(p.adjust(p = simes_pu,
                                            method = "BH") <= alpha)), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
      }
      if(method == "AdaFilter"){
        obj <- adaFilter(p.matrix = p.mat, r = u, type.I.err = "FDR", alpha = alpha/sum(1/(1:m)))
        data.row <- c(length(which(obj$decision == 1)), u, method, alpha)
        repRejection_df <- rbind.data.frame(repRejection_df, data.row)
      }
    }
    
  }
}

colnames(repRejection_df) <- c("rejections", "u", "method", "alpha")
repRejection_df$method <- factor(x = repRejection_df$method, levels = methods)
levels(repRejection_df$method) <- c("ParFilter","Max-ParFilter",
                                    "Adaptive-ParFilter","Adaptive-Max-ParFilter",
                                    "BH-Simes",
                                    "BH-Stouffer", "BH-Fisher", "AdaFilter")

repRejection_df$u <- as.factor(x = repRejection_df$u)
levels(repRejection_df$u) <- c("u = 2","u = 3","u = 4")
repRejection_df$rejections <- as.numeric(repRejection_df$rejections)

repRejection_df$rejections <- log10(repRejection_df$rejections)

repDMDplot <- ggplot(repRejection_df, aes(x=alpha,
                                          y = rejections,
                                          group = method,
                                          linetype=method,
                                          shape=method,
                                          color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,100)[c(1,2,3,4,6,7,8,9)]) +
  scale_color_manual(values = (10:(10+100))[c(1,2,3,4,6,7,8,9)] ) + 
  facet_grid(. ~ u) +
  ylab( TeX(r'($log_{10}(Rejections)$)') ) +
  xlab( TeX(r'($q_{rep}$)') ) +
  theme_bw() +
  theme(legend.title=element_blank())


methods <- c("PF","MaxPF",
             "APF","MaxAPF")
colnames(studyRejectiondf) <- c("rejections", "study", "u", "method", "alpha")
studyRejectiondf$method <- factor(x = studyRejectiondf$method, levels = methods)
levels(studyRejectiondf$method) <- c("ParFilter","Max-ParFilter",
                                     "Adaptive-ParFilter","Adaptive-Max-ParFilter")
studyRejectiondf$u <- as.factor(x = studyRejectiondf$u)
studyRejectiondf$study <- as.factor(studyRejectiondf$study)
levels(studyRejectiondf$u) <- c("u = 2","u = 3","u = 4")
levels(studyRejectiondf$study) <- c("GDS214","GDS563","GDS1956","GDS3027")
studyRejectiondf$rejections <- as.numeric(studyRejectiondf$rejections)
studyRejectiondf$alpha <- as.numeric(studyRejectiondf$alpha)/5
studyDMDplot <- ggplot(studyRejectiondf[studyRejectiondf$u == 'u = 4',], aes(x=alpha,
                                                                             y = rejections,
                                                                             group = method,
                                                                             linetype=method,
                                                                             shape=method,
                                                                             color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))[c(1,2,3,4,6,9)]) +
  scale_color_manual(values = (10:(10+length(methods)))[c(1,2,3,4,6,9)] ) + 
  facet_grid(. ~ study) +
  ylab( TeX(r'(Rejections)') ) +
  xlab( TeX(r'($q_{Study}$)') ) +
  theme_bw() +
  theme(legend.title=element_blank())

pdf(file = "DMDResults44Manual.pdf",  
    width = 10,
    height = 6) 
ggarrange(repDMDplot,
          studyDMDplot,
          ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

# Plotting the histogram
DMD_p1_plot <- ggplot(data.frame(GDS214 = DMD.pvalues[,1]), aes(x=GDS214)) + 
  geom_histogram(aes(y = ..density..)) + geom_vline(xintercept = 0.05, 
                                                    color = "red", size=1)
DMD_p2_plot <- ggplot(data.frame(GDS563 = DMD.pvalues[,2]), aes(x=GDS563)) + 
  geom_histogram(aes(y = ..density..)) + geom_vline(xintercept = 0.05, 
                                                    color = "red", size=1)
DMD_p3_plot <- ggplot(data.frame(GDS1956 = DMD.pvalues[,3]), aes(x=GDS1956)) + 
  geom_histogram(aes(y = ..density..)) + geom_vline(xintercept = 0.05, 
                                                    color = "red", size=1)
DMD_p4_plot <- ggplot(data.frame(GDS3027 = DMD.pvalues[,4]), aes(x=GDS3027)) + 
  geom_histogram(aes(y = ..density..)) + geom_vline(xintercept = 0.05, 
                                                    color = "red", size=1)

pdf(file = "DMDHistogram.pdf",  
    width = 10, 
    height = 3) 
ggarrange(DMD_p1_plot,
          DMD_p2_plot,
          DMD_p3_plot,
          DMD_p4_plot,
          ncol=4, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()
