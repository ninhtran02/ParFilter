# Load the data
load("AIREKO_data.Rdata")
library(CAMT)
library(IHW)
library(adaptMT)
library(adaFilter)
library(ParFilter)

methods <- c("ParFilter", "BH", "Inflated-AdaFilter-BH",
             "AdaFilter-BH", "CAMT", "AdaPT", "IHW", "Naive ParFilter")

q <- 0.05
u <- 3
K <- 3
R_05_list <- rep(list(0), length(methods))
names(R_05_list) <- methods

for(method in methods){
  
  if(method == "ParFilter"){
    R_set <- ParFilter_FDR(p_mat = p.mat, X_list = X_list, u = u, q = q, K = K,
                           method = "Stouffer", adaptive = TRUE, cross_weights = FALSE,
                           lambdas = rep(0.50,K))
    R_05_list[[method]] <- R_set
  }
  
  if(method == "Naive ParFilter"){
    R_set <- ParFilter_FDR(p_mat = p.mat, X_list = X_list, u = u, q = q, K = K,
                           method = "Stouffer", adaptive = TRUE, cross_weights = "naive",
                           lambdas = rep(0.50,K))
    R_05_list[[method]] <- R_set
  }
  
  if(method == "Inflated AdaFilter-BH"){
    R_set <- which(adaFilter(p.matrix = p.mat, r = u,
                             alpha = q/sum(1/(1:m)))$decision == 1)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "AdaFilter-BH"){
    R_set <- which(adaFilter(p.matrix = p.mat, r = u, alpha = q)$decision == 1)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "Oracle"){
    R_set <- oracle_procedure(meta_study = meta_study, u = u, alpha = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "BH"){
    R_set <- BH_procedure(p_mat = p.mat, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "CAMT"){
    R_set <- c()
    R_set <- CAMT_procedure(p_mat = p.mat, u = u, q = q, X_mat = X_spline)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "AdaPT"){
    R_set <- AdaPT_procedure(p_mat = p.mat, X_mat = X_spline, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "IHW"){
    R_set <- IHW_procedure(p_mat = p.mat, X_mat = rowMeans(X), u = u, q = q)
    R_05_list[[method]] <- R_set
  }
}


# Print Genes that the ParFilter can find but not the other methods
other_methods_union_R_set <- unique(unlist(R_05_list[setdiff(methods, "ParFilter")]))
ParFilter_R_set <- R_05_list[["ParFilter"]]
diff_R_set <- setdiff(ParFilter_R_set,other_methods_union_R_set)
diff_p.mat <- p.mat[diff_R_set,]
diff_PC_p <- apply(X = diff_p.mat, MARGIN = 1, FUN = stouffer_fun, u = 3)
sort(diff_PC_p)
names(sort(diff_PC_p))

# Print Genes that the ParFilter did not find compared to other methods
rownames(p_mat)[setdiff(R_05_list$`AdaFilter-BH`,R_05_list$ParFilter)]

# Make the plots
R_05_plot <- ggplot(R_05_dat, aes(x=Method, y= Rejections, fill=Method )) +  
  geom_bar(stat = "identity" ) +
  theme_bw() +
  scale_fill_manual(values = color_ref ) +
  theme(legend.position = "none", axis.text.x = element_text(color = "Black",
                                                             size = 8, angle = 45, vjust = 0.5)) +
  geom_text(aes(label=Rejections), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)

pdf(file = "AIRERejectionsWithNaive05.pdf",
    width = 6, height = 5)
R_05_plot
dev.off()


