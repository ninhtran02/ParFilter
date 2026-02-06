library(ggplot2)
library(latex2exp)
library(ggpubr)
pi1pi1_options <- c(0.75, 0.85, 0.95)
mu_options <- c(0.74, 0.76, 0.78, 0.80, 0.82)
u_n_options <- list(c(2,2),c(2,3),c(3,3),c(3,4),c(4,4),c(3,5),c(4,5),c(5,5))
paral_options <- 1:50
methods <- c("Non-adaptive-ParFilter", "Adaptive-BH",
             "Inflated-ParFilter", "BY", "ParFilter", "BH", "Inflated-AdaFilter-BH",
             "AdaFilter-BH", "CAMT", "AdaPT", "IHW", "CoFilter-BH", "Adaptive-CoFilter-BH",
             "Naive-ParFilter")


FDR_dat <- data.frame()
TPR_dat <- data.frame()
for(pi1pi1 in pi1pi1_options){
  for(mu in mu_options){
    for(u_n in u_n_options){
      paral <- 1
      u <- u_n[1]
      n <- u_n[2]
      file.name <- paste("Independence/uleqn/uleqn m = 5000 0.05/","mu",mu,"setup",u,n,"pi1pi1",pi1pi1,"paral",paral,".RD",sep = "")
      load(file = file.name)
      fdr_dat <- cbind.data.frame(reshape2::melt(FDR_list),pi1pi1,mu,paste(u,n, sep = ""))
      tpr_dat <- cbind.data.frame(reshape2::melt(TPR_list),pi1pi1,mu,paste(u,n, sep = ""))
      fdr_dat$value <- 0
      tpr_dat$value <- 0
      for(paral in paral_options){
        u <- u_n[1]
        n <- u_n[2]
        file.name <- paste("Independence/uleqn/uleqn m = 5000 0.05/","mu",mu,"setup",u,n,"pi1pi1",pi1pi1,"paral",paral,".RD",sep = "")
        load(file = file.name)
        fdr_dat$value <- fdr_dat$value + cbind.data.frame(reshape2::melt(FDR_list),pi1pi1,mu,paste(u,n, sep = ""))$value/max(paral_options)
        tpr_dat$value <- tpr_dat$value + cbind.data.frame(reshape2::melt(TPR_list),pi1pi1,mu,paste(u,n, sep = ""))$value/max(paral_options)
      }
      FDR_dat <- rbind.data.frame(FDR_dat,fdr_dat)
      TPR_dat <- rbind.data.frame(TPR_dat,tpr_dat)
    }
  }
}

names(FDR_dat) <- c("FDR", "Methods", "pi1pi1", "mu", "un")
names(TPR_dat) <- c("TPR", "Methods", "pi1pi1", "mu", "un")

q <- 0.05

FDR_dat$pi1pi1 <- as.factor(x = FDR_dat$pi1pi1)
levels(FDR_dat$pi1pi1) <- c("Least-informative~(delta == 0.75)",
                           "Mildly~informative~(delta == 0.85)",
                           "Most~informative~(delta == 0.95)")
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
                             "Adaptive-CoFilter-BH",
                             "BH", "BY", "CAMT",
                             "CoFilter-BH",
                             "IHW", "Inflated-AdaFilter-BH",
                             "Inflated-ParFilter",
                             "No-Covar-ParFilter",
                             "Non-adaptive-ParFilter",
                             #"Oracle",
                             "ParFilter")

FDR_dat$Methods <- factor(FDR_dat$Methods, levels = c("ParFilter", "BH", "BY", "No-Covar-ParFilter", "AdaPT", "Adaptive-CoFilter-BH", "AdaFilter-BH",  "CAMT", "Inflated-ParFilter",
                                                      "Inflated-AdaFilter-BH", "IHW", "Non-adaptive-ParFilter", "CoFilter-BH",  "Adaptive-BH" #, "Oracle"
))

TPR_dat$pi1pi1 <- as.factor(x = TPR_dat$pi1pi1)
levels(TPR_dat$pi1pi1) <- c("Least-informative~(delta == 0.75)",
                          "Mildly~informative~(delta == 0.85)",
                          "Most~informative~(delta == 0.95)")
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
levels(TPR_dat$Methods) <- c("AdaFilter-BH", "AdaPT", "Adaptive-BH",
                             "Adaptive-CoFilter-BH",
                             "BH", "BY", "CAMT",
                             "CoFilter-BH",
                             "IHW", "Inflated-AdaFilter-BH",
                             "Inflated-ParFilter",
                             "No-Covar-ParFilter",
                             "Non-adaptive-ParFilter",
                             #"Oracle",
                             "ParFilter")
TPR_dat$Methods <- factor(TPR_dat$Methods, levels = c("ParFilter", "BH", "BY", "No-Covar-ParFilter", "AdaPT", "Adaptive-CoFilter-BH", "AdaFilter-BH",  "CAMT", "Inflated-ParFilter",
                                                      "Inflated-AdaFilter-BH", "IHW", "Non-adaptive-ParFilter", "CoFilter-BH",  "Adaptive-BH" #,"Oracle"
))




color_ref <- c("#0006b1", 15, "#795548",  "#ae0000", 11, "#e6d7ff", "#3b444b", 14, "#673AB7", 12, 13, "#808000", "#DE7E5D", "#FF5722"
               #,16
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

shape_ref <- c(7, 2, 9, 14, 1, 12, 0, 3, 10, 5, 4, 11, 13, 8
               #, 6
)
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

alpha_ref <- c(1, 0.75, 0.75, 1, 0.75, 0.75, 0.9, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75
               #, 0.75
)
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

line_ref <- c(1, 11, 10, 1, 14, 12, 3, 9, 5, 6, 7, 4, 8, 13#, 2
)
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

FDR_plot <- ggplot(data=FDR_dat[FDR_dat$un %in% c("2 / '['*2*']'",
                                                  "3 / '['*3*']'",
                                                  "4 / '['*4*']'",
                                                  "5 / '['*5*']'"),], aes(x=mu,
                                                                          y = FDR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods,
                                                                          alpha=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=shape_ref) +
  scale_color_manual(values = color_ref ) +
  scale_linetype_manual(values = line_ref) +
  scale_alpha_manual(values = alpha_ref) +
  facet_grid(un ~ pi1pi1, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  geom_hline(yintercept = q, linetype='dotted', size = 1) +
  guides(color = guide_legend(direction = "horizontal")) +
  theme(legend.title=element_blank(), legend.box = "horizontal")

TPR_plot <- ggplot(data=TPR_dat[TPR_dat$un %in% c("2 / '['*2*']'",
                                                  "3 / '['*3*']'",
                                                  "4 / '['*4*']'",
                                                  "5 / '['*5*']'"),], aes(x=mu,
                                                                          y = TPR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods,
                                                                          alpha=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=shape_ref) +
  scale_color_manual(values = color_ref ) +
  scale_linetype_manual(values = line_ref) +
  scale_alpha_manual(values = alpha_ref) +
  facet_grid(un ~ pi1pi1, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank(), legend.box = "horizontal")

#pdf(file = "u = n Results.pdf",
#    width = 9.6, height = 9.6)#12.8)
pdf(file = "u_equal_n_Results_uleqn.pdf",
    width = 9.6, height = 13.6)
ggarrange(FDR_plot,TPR_plot,ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()

FDR_plot <- ggplot(data=FDR_dat[FDR_dat$un %in% c("2 / '['*3*']'",
                                                  "3 / '['*4*']'",
                                                  "3 / '['*5*']'",
                                                  "4 / '['*5*']'"),], aes(x=mu,
                                                                          y = FDR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods,
                                                                          alpha=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=shape_ref) +
  scale_color_manual(values = color_ref ) +
  scale_linetype_manual(values = line_ref) +
  scale_alpha_manual(values = alpha_ref) +
  facet_grid(un ~ pi1pi1, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  guides(color = guide_legend(direction = "horizontal")) +
  theme(legend.title=element_blank(), legend.position = "bottom")

TPR_plot <- ggplot(data=TPR_dat[TPR_dat$un %in% c("2 / '['*3*']'",
                                                  "3 / '['*4*']'",
                                                  "3 / '['*5*']'",
                                                  "4 / '['*5*']'"),], aes(x=mu,
                                                                          y = TPR,
                                                                          group = Methods,
                                                                          linetype=Methods,
                                                                          shape=Methods,
                                                                          color=Methods,
                                                                          alpha=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=shape_ref) +
  scale_color_manual(values = color_ref ) +
  scale_linetype_manual(values = line_ref) +
  scale_alpha_manual(values = alpha_ref) +
  facet_grid(un ~ pi1pi1, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($\xi$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  guides(color = guide_legend(direction = "horizontal")) +
  theme(legend.title=element_blank(), legend.position = "bottom")

pdf(file = "u_less_n_Results_uleqn.pdf",
    width = 9.6, height = 13.6)
ggarrange(FDR_plot,TPR_plot,ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()
