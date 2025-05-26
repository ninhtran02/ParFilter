setwd("~/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Paper Simulations/Covariate-Assisted")
library(ggplot2)
library(latex2exp)
library(ggpubr)
xcoef_options <- c(0, 1, 1.5)
mu_options <- c(0.74, 0.76, 0.78, 0.80, 0.82)
u_n_options <- list(c(2,2),c(2,3),c(3,3),c(3,4),c(4,4),c(3,5),c(4,5),c(5,5))
methods <- c("ParFilter", "BH", "Inflated-AdaFilter-BH",
             "AdaFilter-BH", "CAMT", "AdaPT", "IHW", 
             "Non-adaptive-ParFilter", "Adaptive-BH",
             "Inflated-ParFilter", "BY")


FDR_dat <- data.frame()
TPR_dat <- data.frame()
for(xcoef in xcoef_options){
  for(mu in mu_options){
    for(u_n in u_n_options){
      u <- u_n[1]
      n <- u_n[2]
      file.name <- paste("NegativeDependence/","mu",mu,"setup","xcoef",xcoef,u,n,".RD",sep = "")
      load(file = file.name)
      fdr_dat <- cbind.data.frame(reshape2::melt(FDR_list),xcoef,mu,paste(u,n, sep = ""))
      FDR_dat <- rbind.data.frame(FDR_dat,fdr_dat)
      tpr_dat <- cbind.data.frame(reshape2::melt(TPR_list),xcoef,mu,paste(u,n, sep = ""))
      TPR_dat <- rbind.data.frame(TPR_dat,tpr_dat)
    }
  }
}

names(FDR_dat) <- c("FDR", "Methods", "xcoef", "mu", "un")
names(TPR_dat) <- c("TPR", "Methods", "xcoef", "mu", "un")

FDR_dat$xcoef <- as.factor(x = FDR_dat$xcoef)
levels(FDR_dat$xcoef) <- c("Non-informative~(gamma[1] == 0)",
                           "Mildly~informative~(gamma[1] == 1.0)",
                           "Most~informative~(gamma[1] == 1.5)")
FDR_dat$un <- as.factor(x = FDR_dat$un)
levels(FDR_dat$un) <- c("2/[2]","2/[3]","3/[3]","3/[4]","3/[5]","4/[4]","4/[5]","5/[5]")
levels(FDR_dat$un) <- c("2 / '['*2*']'",
                        "2 / '['*3*']'",
                        "3 / '['*3*']'",
                        "3 / '['*4*']'",
                        "3 / '['*5*']'",
                        "4 / '['*4*']'",
                        "4 / '['*5*']'",
                        "5 / '['*5*']'")

FDR_dat$Methods <- as.factor(x = FDR_dat$Methods)
levels(FDR_dat$Methods) <- c("AdaFilter-BH", "AdaPT", "Adaptive-BH", 
                             "BH", "BY", "CAMT", "IHW", "Inflated-AdaFilter-BH",
                             "Inflated-ParFilter", "Non-adaptive-ParFilter", "ParFilter")

TPR_dat$xcoef <- as.factor(x = TPR_dat$xcoef)
levels(TPR_dat$xcoef) <- c("Non-informative~(gamma[1] == 0)",
                           "Mildly~informative~(gamma[1] == 1.0)",
                           "Most~informative~(gamma[1] == 1.5)")
TPR_dat$un <- as.factor(x = TPR_dat$un)
levels(TPR_dat$un) <- c("2/[2]","2/[3]","3/[3]","3/[4]","3/[5]","4/[4]","4/[5]","5/[5]")
levels(TPR_dat$un) <- c("2 / '['*2*']'",
                        "2 / '['*3*']'",
                        "3 / '['*3*']'",
                        "3 / '['*4*']'",
                        "3 / '['*5*']'",
                        "4 / '['*4*']'",
                        "4 / '['*5*']'",
                        "5 / '['*5*']'")

TPR_dat$Methods <- as.factor(x = TPR_dat$Methods)
levels(TPR_dat$Methods) <-  c("AdaFilter-BH", "AdaPT", "Adaptive-BH", 
                              "BH", "BY", "CAMT", "IHW", "Inflated-AdaFilter-BH",
                              "Inflated-ParFilter", "Non-adaptive-ParFilter", "ParFilter")

#  c("AdaFilter-BH", "AdaPT", "BY", "CAMT",
#                             "IHW", "Inflated-AdaFilter-BH", "Inflated-ParFilter", "ParFilter")
color_ref <- c(17, 11, # "AdaFilter-BH", "AdaPT"
               2, # "Adaptive-BH"
               15, # "BH"
               3, # "BY"
               14, 13, 12, #"CAMT", "IHW", "Inflated-AdaFilter-BH"
               16, #"Inflated-ParFilter"
               4, # "Non-adaptive-ParFilter"
               10) # "ParFilter"

shape_ref <- c(0, 1, # "AdaFilter-BH", "AdaPT"
               8, # "Adaptive-BH"
               2,  # "BH"
               9, # "BY"
               3,4,5, #"CAMT", "IHW", "Inflated-AdaFilter-BH"
               10, #"Inflated-ParFilter"
               11, # "Non-adaptive-ParFilter"
               #6, # "Oracle"
               7) # "ParFilter"

FDR_plot <- ggplot(data=FDR_dat[FDR_dat$un %in% c("2 / '['*2*']'",
                                                  "3 / '['*3*']'",
                                                  "4 / '['*4*']'",
                                                  "5 / '['*5*']'"),], aes(x=mu,
                                                                          y = FDR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=shape_ref) +
  scale_color_manual(values = color_ref ) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

TPR_plot <- ggplot(data=TPR_dat[TPR_dat$un %in% c("2 / '['*2*']'",
                                                  "3 / '['*3*']'",
                                                  "4 / '['*4*']'",
                                                  "5 / '['*5*']'"),], aes(x=mu,
                                                                          y = TPR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=shape_ref) +
  scale_color_manual(values = color_ref ) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

pdf(file = "dep_u_equal_n_Results_full.pdf",
    width = 9.6, height = 9.6)#12.8)
ggarrange(FDR_plot,TPR_plot,ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()

#FDR_plot <- ggplot(data=FDR_dat[FDR_dat$un %in% c("2 / '['*3*']'",
#                                                  "3 / '['*4*']'",
#                                                  "3 / '['*5*']'",
#                                                  "4 / '['*5*']'"),], aes(x=mu,
#                                                                          y = FDR,
#                                                                          group = Methods,
#                                                                          linetype=Methods,
#                                                                          shape=Methods,
#                                                                          color=Methods)) +
#  geom_line() +
#  geom_point() +
#  scale_shape_manual(values=shape_ref) +
#  scale_color_manual(values = 10:(10+length(methods)) ) +
#  facet_grid(un ~ xcoef, labeller = label_parsed) +
#  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
#  ylab(TeX(r'(FDR$_{rep}$)')) +
#  xlab( TeX(r'($\mu$)') ) +
#  theme_bw() +
#  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
#  theme(legend.title=element_blank())
#
#TPR_plot <- ggplot(data=TPR_dat[TPR_dat$un %in% c("2 / '['*3*']'",
#                                                  "3 / '['*4*']'",
#                                                  "3 / '['*5*']'",
#                                                  "4 / '['*5*']'"),], aes(x=mu,
#                                                                          y = TPR,
#                                                                          group = Methods,
#                                                                          linetype=Methods,
#                                                                          shape=Methods,
#                                                                          color=Methods)) +
#  geom_line() +
#  geom_point() +
#  scale_shape_manual(values=shape_ref) +
#  scale_color_manual(values = 10:(10+length(methods)) ) +
#  facet_grid(un ~ xcoef, labeller = label_parsed) +
#  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
#  ylab(TeX(r'(TPR$_{rep}$)')) +
#  xlab( TeX(r'($\mu$)') ) +
#  theme_bw() +
#  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
#  theme(legend.title=element_blank())
#
#pdf(file = "dep u < n Results.pdf 2025",
#    width = 9.6, height = 9.6)
#ggarrange(FDR_plot,TPR_plot,ncol = 1, nrow = 2,
#          common.legend = TRUE, legend = "bottom")
#dev.off()
#
#
#