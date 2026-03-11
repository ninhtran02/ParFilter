setwd("~/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Paper Simulations/Covariate-Assisted/Second Round Codes")
library(ggplot2)
library(latex2exp)
library(ggpubr)
xcoef_options <- c(0.0, 1.0, 1.25)
mu_options <- c(0.74, 0.76, 0.78, 0.80, 0.82)
u_n_options <- list(c(5,5),c(6,6),c(7,7),c(8,8))
paral_options <- 1:50
methods <- c("ParFilter_0", "ParFilter_1",
             "ParFilter_2", "ParFilter_3",
             "No_Covar_ParFilter_0", "No_Covar_ParFilter_1",
             "No_Covar_ParFilter_2", "No_Covar_ParFilter_3")


FDR_dat <- data.frame()
TPR_dat <- data.frame()
for(xcoef in xcoef_options){
  for(mu in mu_options){
    for(u_n in u_n_options){
      paral <- 1
      u <- u_n[1]
      n <- u_n[2]
      file.name <- paste("Independence/PF/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
      load(file = file.name)
      fdr_dat <- cbind.data.frame(reshape2::melt(FDR_list),xcoef,mu,paste(u,n, sep = ""))
      tpr_dat <- cbind.data.frame(reshape2::melt(TPR_list),xcoef,mu,paste(u,n, sep = ""))
      fdr_dat$value <- 0
      tpr_dat$value <- 0
      for(paral in paral_options){
        u <- u_n[1]
        n <- u_n[2]
        file.name <- paste("Independence/PF/","mu",mu,"setup",u,n,"xcoef",xcoef,"paral",paral,".RD",sep = "")
        load(file = file.name)
        fdr_dat$value <- fdr_dat$value + cbind.data.frame(reshape2::melt(FDR_list),xcoef,mu,paste(u,n, sep = ""))$value/max(paral_options)
        tpr_dat$value <- tpr_dat$value + cbind.data.frame(reshape2::melt(TPR_list),xcoef,mu,paste(u,n, sep = ""))$value/max(paral_options)
      }
      FDR_dat <- rbind.data.frame(FDR_dat,fdr_dat)
      TPR_dat <- rbind.data.frame(TPR_dat,tpr_dat)
    }
  }
}

names(FDR_dat) <- c("FDR", "Methods", "xcoef", "mu", "un")
names(TPR_dat) <- c("TPR", "Methods", "xcoef", "mu", "un")

q <- 0.05

FDR_dat$xcoef <- as.factor(x = FDR_dat$xcoef)
levels(FDR_dat$xcoef) <- c("Non-informative~(gamma[1] == 0)",
                           "Mildly~informative~(gamma[1] == 1.0)",
                           "Most~informative~(gamma[1] == 1.25)")
FDR_dat$un <- as.factor(x = FDR_dat$un)
levels(FDR_dat$un) <- c("5/[5]","6/[6]","7/[7]","8/[8]")
levels(FDR_dat$un) <- c("5 / '['*5*']'",
                        "6 / '['*6*']'",
                        "7 / '['*7*']'",
                        "8 / '['*8*']'")

FDR_dat$Methods <- factor(x = FDR_dat$Methods,
                            levels = c("ParFilter_0",
                                       "No_Covar_ParFilter_0",
                                       "ParFilter_1",
                                       "No_Covar_ParFilter_1",
                                       "ParFilter_2",
                                       "No_Covar_ParFilter_2",
                                       "ParFilter_3",
                                       "No_Covar_ParFilter_3"))
levels(FDR_dat$Methods) <- c("ParFilter[n]",
                             "No-Covar-ParFilter[n]",
                             "ParFilter[n-1]",
                             "No-Covar-ParFilter[n-1]",
                             "ParFilter[n-2]",
                             "No-Covar-ParFilter[n-2]",
                             "ParFilter[n-3]",
                             "No-Covar-ParFilter[n-3]")

#FDR_dat$Methods <- factor(FDR_dat$Methods, levels = c("ParFilter", "BH", "BY", "No-Covar-ParFilter", "AdaPT", "Adaptive-CoFilter-BH", "AdaFilter-BH",  "CAMT", "Inflated-ParFilter",
#                                                      "Inflated-AdaFilter-BH", "IHW", "Non-adaptive-ParFilter", "CoFilter-BH",  "Adaptive-BH" #, "Oracle"
#))

TPR_dat$xcoef <- as.factor(x = TPR_dat$xcoef)
levels(TPR_dat$xcoef) <- c("Non-informative~(gamma[1] == 0)",
                           "Mildly~informative~(gamma[1] == 1.0)",
                           "Most~informative~(gamma[1] == 1.25)")
TPR_dat$un <- as.factor(x = TPR_dat$un)
levels(TPR_dat$un) <- c("5/[5]","6/[6]","7/[7]","8/[8]")
levels(TPR_dat$un) <- c("5 / '['*5*']'",
                        "6 / '['*6*']'",
                        "7 / '['*7*']'",
                        "8 / '['*8*']'")

TPR_dat$Methods <- factor(x = TPR_dat$Methods,
                          levels = c("ParFilter_0",
                                     "No_Covar_ParFilter_0",
                                     "ParFilter_1",
                                     "No_Covar_ParFilter_1",
                                     "ParFilter_2",
                                     "No_Covar_ParFilter_2",
                                     "ParFilter_3",
                                     "No_Covar_ParFilter_3"))
levels(TPR_dat$Methods) <- c("ParFilter[n]",
                             "No-Covar-ParFilter[n]",
                             "ParFilter[n-1]",
                             "No-Covar-ParFilter[n-1]",
                             "ParFilter[n-2]",
                             "No-Covar-ParFilter[n-2]",
                             "ParFilter[n-3]",
                             "No-Covar-ParFilter[n-3]")
#levels(TPR_dat$Methods) <- c("AdaFilter-BH", "AdaPT", "Adaptive-BH", 
#                             "Adaptive-CoFilter-BH",
#                             "BH", "BY", "CAMT", 
#                             "CoFilter-BH",
#                             "IHW", "Inflated-AdaFilter-BH",
#                             "Inflated-ParFilter",
#                             "No-Covar-ParFilter",
#                             "Non-adaptive-ParFilter", 
#                             #"Oracle", 
#                             "ParFilter")
#TPR_dat$Methods <- factor(TPR_dat$Methods, levels = c("ParFilter", "BH", "BY", "No-Covar-ParFilter", "AdaPT", "Adaptive-CoFilter-BH", "AdaFilter-BH",  "CAMT", "Inflated-ParFilter",
#                                                      "Inflated-AdaFilter-BH", "IHW", "Non-adaptive-ParFilter", "CoFilter-BH",  "Adaptive-BH" #,"Oracle"
#))




color_ref <- c("#0006b1", "#ae0000",
               "#00A676",  # teal-green
               "#FF8C00",  # dark orange
               "#7B2CBF",  # purple
               "#2D6A4F",  # forest green
               "#C9A227",  # golden yellow
               "#00B4D8"  # sky cyan
               )

#color_ref <- c("#3b444b", 11, # "AdaFilter-BH", "AdaPT"
#               "#FF5722", # "Adaptive-BH"
#               "#e6d7ff", #Adaptive-CoFilter-BH
#               15, # "BH"
#               "#795548", # "BY"
#               14, #"CAMT",
#               "#DE7E5D", # "CoFilter-BH"
#               13, 12, # "IHW", #"Inflated-AdaFilter-BH"
#               "#673AB7", #"Inflated-ParFilter"
#               "#ae0000", # "Naive-ParFilter"
#               "#808000", # "Non-adaptive-ParFilter"
#               16, # "Oracle"
#               "#0006b1") # "ParFilter"

shape_ref <- c(7, 14, 2, 9, 1, 12, 0, 3)
#shape_ref <- c(0, 1, # "AdaFilter-BH", "AdaPT"
#               8, # "Adaptive-BH"
#               12, # #Adaptive-CoFilter-BH
#               2,  # "BH"
#               9, # "BY"
#               3, #"CAMT"
#               13, # CoFilter-BH
#               4,5, # "IHW", "Inflated-AdaFilter-BH"
#               10, #"Inflated-ParFilter"
#               14,  # "No-Covar-ParFilter"
#               11, # "Non-adaptive-ParFilter"
#               6, # "Oracle"
#               7) # "ParFilter"

alpha_ref <- c(1, 1, 1, 1, 1, 1, 1, 1)
#alpha_ref <- c(0.90, 0.75, # "AdaFilter-BH", "AdaPT"
#               0.75, # "Adaptive-BH"
#               0.75, # Adaptive-CoFilter-BH
#               0.75, # "BH"
#               0.75, # "BY"
#               0.75, #"CAMT"
#               0.75, # CoFilter-BH
#               0.75, 0.75, # "IHW", "Inflated-AdaFilter-BH"
#               0.75, #"Inflated-ParFilter"
#               1, # "No-Covar-ParFilter"
#               0.75, # "Non-adaptive-ParFilter"
#               0.75, # "Oracle"
#               1) # "ParFilter"

line_ref <- c(1, 1, 11, 10, 1, 14, 12, 3)
#line_ref <- c(3, # "AdaFilter-BH"
#              14, # "AdaPT"
#              13, # "Adaptive-BH"
#              12, # Adaptive-CoFilter-BH
#              11, # BH
#              10, #BY
#              9,  # CAMT
#              8,  # CoFilter-BH
#              7, # "IHW"
#              6,  #"Inflated-AdaFilter-BH"
#              5,  #"Inflated-ParFilter"
#              1,  # "No-Covar-ParFilter"
#              4,  # "Non-adaptive-ParFilter"
#              2,  # "Oracle"
#              1) # "ParFilter"

FDR_plot <- ggplot(data=FDR_dat[FDR_dat$un %in% c("5 / '['*5*']'",
                                                  "6 / '['*6*']'", 
                                                  "7 / '['*7*']'",
                                                  "8 / '['*8*']'"),], aes(x=mu,
                                                                          y = FDR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(breaks = brks, values = shape_ref, labels = labs) + 
  scale_color_manual(breaks = brks, values = color_ref, labels = labs) +
  scale_linetype_manual(breaks = brks, values = line_ref, labels = labs) +
  #scale_alpha_manual(values = alpha_ref) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  geom_hline(yintercept = q, linetype='dotted', size = 1) +
  guides(color = guide_legend(direction = "horizontal")) + 
  theme(legend.title=element_blank(), legend.box = "horizontal")

TPR_plot <- ggplot(data=TPR_dat[TPR_dat$un %in% c("5 / '['*5*']'",
                                                  "6 / '['*6*']'", 
                                                  "7 / '['*7*']'",
                                                  "8 / '['*8*']'"),], aes(x=mu,
                                                                          y = TPR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(breaks = brks, values = shape_ref, labels = labs) + 
  scale_color_manual(breaks = brks, values = color_ref, labels = labs) +
  scale_linetype_manual(breaks = brks, values = line_ref, labels = labs) +
  #scale_alpha_manual(values = alpha_ref) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank(), legend.box = "horizontal") 

pdf(file = "PF_comparison_large_ALL.pdf",
    width = 9.6, height = 13.6)
ggarrange(FDR_plot,TPR_plot,ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()



# Refinding the data
FDR_dat_subset <- FDR_dat
FDR_dat_subset <- FDR_dat_subset[FDR_dat_subset$xcoef == "Most~informative~(gamma[1] == 1.25)",]
TPR_dat_subset <- TPR_dat
TPR_dat_subset <- TPR_dat_subset[TPR_dat_subset$xcoef == "Most~informative~(gamma[1] == 1.25)",]
#
labs <- c(expression(ParFilter[n]),
          expression(No-Covar-ParFilter[n]),
          expression(ParFilter[n-1]),
          expression(No-Covar-ParFilter[n-1]),
          expression(ParFilter[n-2]),
          expression(No-Covar-ParFilter[n-2]),
          expression(ParFilter[n-3]),
          expression(No-Covar-ParFilter[n-3]))
brks <- levels(FDR_dat$Methods)
#
FDR_plot <- ggplot(data=FDR_dat_subset, aes(x=mu,
                                            y = FDR,
                                            group = Methods,
                                            linetype=Methods,
                                            shape=Methods,
                                            color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(breaks = brks, values = shape_ref, labels = labs) + 
  scale_color_manual(breaks = brks, values = color_ref, labels = labs) +
  scale_linetype_manual(breaks = brks, values = line_ref, labels = labs) +
  #scale_alpha_manual(values = alpha_ref) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  geom_hline(yintercept = q, linetype='dotted', size = 1) +
  guides(color = guide_legend(direction = "horizontal")) + 
  theme(legend.title=element_blank(), legend.box = "horizontal")

TPR_plot <- ggplot(data=TPR_dat_subset, aes(x=mu,
                                            y = TPR,
                                            group = Methods,
                                            linetype=Methods,
                                            shape=Methods,
                                            color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(breaks = brks, values = shape_ref, labels = labs) + 
  scale_color_manual(breaks = brks, values = color_ref, labels = labs) +
  scale_linetype_manual(breaks = brks, values = line_ref, labels = labs) +
  #scale_alpha_manual(values = alpha_ref) +
  facet_grid(xcoef ~ un, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank(), legend.box = "horizontal") 

pdf(file = "PF_comparison_large.pdf",
    width = 9.6, height = 3.2)
ggarrange(TPR_plot,ncol = 1, nrow = 1,
          common.legend = TRUE, legend = "bottom")
dev.off()
