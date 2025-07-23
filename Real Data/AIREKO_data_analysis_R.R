# Load the data
#load("AIREKO_data.Rdata")

#gene.names <- rownames(p.mat)
#save(gene.names, file = "genenames.Rdata")  

library(CAMT)
library(IHW)
library(adaptMT)
library(adaFilter)
library(ParFilter)

methods <- c("ParFilter", "BH", "Inflated-AdaFilter-BH",
             "AdaFilter-BH", "CAMT", "AdaPT", "IHW",
             "CoFilter-BH", "No-Covar-ParFilter")

q <- 0.05
u <- 3
K <- 3
R_05_list <- rep(list(0), length(methods))
names(R_05_list) <- methods

m <- nrow(p.mat)

for(method in methods){
  
  if(method == "ParFilter"){
    R_set <- ParFilter_FDR(p_mat = p.mat, X_list = X_list, u = u, q = q, K = K,
                           method = "Stouffer", adaptive = TRUE, cross_weights = FALSE,
                           lambdas = rep(0.50,K))
    R_05_list[[method]] <- R_set
  }
  
  if(method == "No-Covar-ParFilter"){
    R_set <- ParFilter_FDR(p_mat = p.mat, X_list = X_list, u = u, q = q, K = K,
                           method = "Stouffer", adaptive = TRUE, cross_weights = "naive",
                           lambdas = rep(0.50,K))
    R_05_list[[method]] <- R_set
  }
  
  if(method == "Inflated-AdaFilter-BH"){
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
  
  if(method == "CoFilter-BH"){ 
    R_set <- CoFilter_procedure(p_mat = p.mat, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "Adaptive-CoFilter-BH"){ 
    R_set <- Adaptive_CoFilter_procedure(p_mat = p.mat, u = u, q = q)
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
    R_set <- IHW_procedure(p_mat = p.mat, X_mat = rowMeans(X_spline), u = u, q = q)
    R_05_list[[method]] <- R_set
  }
}

#save(R_05_list, file = "AIREKO_results.Rdata")
#load("AIREKO_results.Rdata")

# Print Genes that the ParFilter can find but not the other methods
other_methods_union_R_set <- unique(unlist(R_05_list[setdiff(methods, "ParFilter")]))
ParFilter_R_set <- R_05_list[["ParFilter"]]
diff_R_set <- setdiff(ParFilter_R_set,other_methods_union_R_set)
length(diff_R_set)
diff_p.mat <- p.mat[diff_R_set,]
diff_PC_p <- apply(X = diff_p.mat, MARGIN = 1, FUN = stouffer_fun, u = 3)
sort(diff_PC_p)
names(sort(diff_PC_p))

# Print Genes that the ParFilter and No-Covar-ParFilter can find but not the other methods
methods_of_interest <- setdiff(setdiff(methods, "ParFilter"),"No-Covar-ParFilter")
other_methods_union_R_set <- unique(unlist(R_05_list[methods_of_interest]))
Both_ParFilter_R_set <- unique(c(R_05_list[["ParFilter"]],R_05_list[["No-Covar-ParFilter"]]))
Both_ParFilter_diff_R_set <- setdiff(Both_ParFilter_R_set ,other_methods_union_R_set)
length(Both_ParFilter_diff_R_set)
Both_ParFilter_diff_p.mat <- p.mat[Both_ParFilter_diff_R_set,]
Both_ParFilter_diff_PC_p <- apply(X = Both_ParFilter_diff_p.mat, MARGIN = 1, FUN = stouffer_fun, u = 3)
sort(Both_ParFilter_diff_PC_p)
names(sort(Both_ParFilter_diff_PC_p))[1:15]

rownames(p.mat)[setdiff(ParFilter_R_set, R_05_list$`No-Covar-ParFilter`)]

rownames(p.mat)[setdiff(R_05_list$`AdaFilter-BH`,ParFilter_R_set)]
rownames(p.mat)[setdiff(R_05_list$BH,ParFilter_R_set)]
rownames(p.mat)[setdiff(R_05_list$`CoFilter-BH`,ParFilter_R_set)]
rownames(p.mat)[setdiff(R_05_list$`Inflated-AdaFilter-BH`,ParFilter_R_set)]
rownames(p.mat)[setdiff(R_05_list$CAMT,ParFilter_R_set)]
rownames(p.mat)[setdiff(R_05_list$AdaPT,ParFilter_R_set)]
rownames(p.mat)[setdiff(R_05_list$IHW,ParFilter_R_set)]

methods_index <- c(1,2,5,7,8,9,10,12,15)
color_ref <- c("#3b444b", 11, # "AdaFilter-BH", "AdaPT"
               "#FF5722", # "Adaptive-BH"
               "#e6d7ff", #Adaptive-CoFilter-BH
               15, # "BH"
               "#795548", # "BY"
               14, #"CAMT",
               "#DE7E5D", # "CoFilter-BH"
               13, 12, # "IHW", #"Inflated-AdaFilter-BH"
               "#673AB7", #"Inflated-ParFilter"
               "#ae0000", # "No-Covar-ParFilter"
               "#808000", # "Non-adaptive-ParFilter"
               16, # "Oracle"
               "#0006b1") # "ParFilter"

alpha_ref <- c(0.70, 0.70, # "AdaFilter-BH", "AdaPT"
               0.70, # "Adaptive-BH"
               0.70, # #Adaptive-CoFilter-BH
               0.70, # "BH"
               0.70, # "BY"
               0.70, #"CAMT"
               0.70, # CoFilter-BH
               0.70, 0.70, # "IHW", "Inflated-AdaFilter-BH"
               0.70, #"Inflated-ParFilter"
               1, # "No-Covar-ParFilter"
               0.70, # "Non-adaptive-ParFilter"
               0.70, # "Oracle"
               1) # "ParFilter"

color_ref <- color_ref[methods_index]
alpha_ref <- alpha_ref[methods_index]

R_05_dat <- data.frame(Method = names(R_05_list), Rejections = unlist(lapply(R_05_list, length)))
R_05_dat$Method <- as.factor(R_05_dat$Method)
# Make the plots
R_05_plot <- ggplot(R_05_dat, aes(x=Method, y= Rejections, fill=Method )) +  
  geom_bar(stat = "identity" ) +
  theme_bw() +
  scale_fill_manual(values = color_ref ) +
  scale_alpha_manual(values = alpha_ref) +
  theme(legend.position = "none", axis.text.x = element_text(color = "Black",
                                                             size = 8, angle = 45, vjust = 0.5)) +
  geom_text(aes(label=Rejections), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)

pdf(file = "AIRERejectionsWithNaiveCoFilter05NoCovarParFilter.pdf",
    width = 6, height = 5)
R_05_plot
dev.off()


