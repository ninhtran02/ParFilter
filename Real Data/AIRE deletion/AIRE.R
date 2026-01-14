setwd("/Users/ninht/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Paper Simulations/Covariate-Assisted/Second Round Codes/RealData2")
# Load the data
load("AIREKO_data.Rdata")
p_mat <- p.mat

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
lambda <- 0.5
R_05_list <- rep(list(0), length(methods))
names(R_05_list) <- methods
PF_covars <- cbind(splines::bs(X[,2], df = 6))

m <- nrow(p_mat)

for(method in methods){
  
  if(method == "ParFilter"){
    R_set <- ParFilter(P = p_mat, u = 3, X_list = list(PF_covars,
                                                       PF_covars,
                                                       PF_covars),
                       K = 3, direction = "negative",
                       q = 0.05, lambda_vec =  rep(lambda,K))
    R_05_list[[method]] <- R_set
  }
  
  if(method == "No-Covar-ParFilter"){
    X_list_no_covar <- list(rep(0,m),rep(0,m),rep(0,m))
    R_set <- ParFilter(P = p_mat, u = 3, X_list = X_list_no_covar,
                       K = 3, direction = "negative",
                       q = 0.05, lambda_vec =  rep(lambda,K))
    R_05_list[[method]] <- R_set
  }
  
  if(method == "Inflated-AdaFilter-BH"){
    R_set <- AdaFilter_procedure(p_mat = p_mat, u = u, q = q/sum(1/(1:m)))
    R_05_list[[method]] <- R_set
  }
  
  if(method == "AdaFilter-BH"){
    R_set <- AdaFilter_procedure(p_mat = p_mat, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "CoFilter-BH"){ 
    R_set <- CoFilter_procedure(p_mat = p_mat, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "Adaptive-CoFilter-BH"){ 
    R_set <- Adaptive_CoFilter_procedure(p_mat = p_mat, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "BH"){
    R_set <- BH_procedure(p_mat = p_mat, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "CAMT"){
    R_set <- c()
    R_set <- CAMT_procedure(p_mat = p_mat, u = u, q = q, X_mat = PF_covars)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "AdaPT"){
    R_set <- AdaPT_procedure(p_mat = p_mat, X_mat = PF_covars, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
  
  if(method == "IHW"){
    IHW_covars <- X[,2]
    R_set <- IHW_procedure(p_mat = p_mat, X_mat = IHW_covars, u = u, q = q)
    R_05_list[[method]] <- R_set
  }
}

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
  #geom_text(aes(label=Rejections), vjust=1.6, color="white",
  #          position = position_dodge(0.9), size=3.5)
  geom_text(aes(label=Rejections), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3.5)

pdf(file = "AIRE_lambda050.pdf",
    width = 6, height = 5)
R_05_plot
dev.off()

reference_PC_p_values <- apply(X = p_mat, MARGIN = 1, FUN = stouffer_fun, u = 3)
names(reference_PC_p_values) <- rownames(p_mat)

not_PF_R_set <- unique(c(R_05_list$BH,
                         R_05_list$`Inflated-AdaFilter-BH`,
                         R_05_list$`AdaFilter-BH`,
                         R_05_list$CAMT,
                         R_05_list$AdaPT,
                         R_05_list$IHW,
                         R_05_list$`CoFilter-BH`,
                         R_05_list$`No-Covar-ParFilter`))
PF_R_set <- unique(c(R_05_list$ParFilter))

length(sort(reference_PC_p_values[setdiff(PF_R_set,not_PF_R_set)]))
# TOP 15 genes by PF
cbind(sort(reference_PC_p_values[setdiff(PF_R_set,not_PF_R_set)])[1:15])


### PF and No Covar
not_PF_R_set <- unique(c(R_05_list$BH,
                         R_05_list$`Inflated-AdaFilter-BH`,
                         R_05_list$`AdaFilter-BH`,
                         R_05_list$CAMT,
                         R_05_list$AdaPT,
                         R_05_list$IHW,
                         R_05_list$`CoFilter-BH`))
PF_R_set <- unique(c(R_05_list$ParFilter,R_05_list$`No-Covar-ParFilter`))

length(sort(reference_PC_p_values[setdiff(PF_R_set,not_PF_R_set)]))
# TOP 15 genes by PF and No-covar-PF
cbind(sort(reference_PC_p_values[setdiff(PF_R_set,not_PF_R_set)])[1:15])




