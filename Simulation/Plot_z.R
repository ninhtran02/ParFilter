setwd("~/Paper Simulations/Covariate-Assisted")
library(ggplot2)
library(latex2exp)
library(ggpubr)
xcoef_options <- c(0, 1.0, 1.5)
mu_options <- c(2.0, 2.2, 2.4, 2.6, 2.8)
u_n_options <- list(c(2,2),c(2,3),c(3,3),c(3,4),c(4,4),c(3,5),c(4,5),c(5,5))
methods <- c("ParFilter", "BH", "Inflated-AdaFilter-BH",
             "AdaFilter-BH", "CAMT", "AdaPT", "IHW", "Oracle", "repfdr")

FDR_dat <- data.frame()
TPR_dat <- data.frame()
for(xcoef in xcoef_options){
  for(mu in mu_options){
    for(u_n in u_n_options){
      u <- u_n[1]
      n <- u_n[2]
      file.name <- paste("Independencez/","mu",mu,"setup",u,n,"xcoef",xcoef,".RD",sep = "")
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
#levels(FDR_dat$Methods) <- c("AdaFilter-BH", "AdaPT", "BH", "CAMT",                 
#                             "IHW", "Inflated-AdaFilter-BH", "Inflated-ParFilter", "ParFilter")

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
#levels(TPR_dat$Methods) <- c("AdaFilter-BH", "AdaPT", "BH", "CAMT",                 
#                             "IHW", "Inflated-AdaFilter-BH", "Inflated-ParFilter", "ParFilter")

color_ref <- c(17, 11, 15, 14, 13, 12, 16, 10, 19)

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
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = color_ref ) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\mu$)') ) +
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
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = color_ref ) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

pdf(file = "z u = n Results.pdf",
    width = 9.6, height = 9.6)
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
                                                                                  color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

TPR_plot <- ggplot(data=TPR_dat[TPR_dat$un %in% c("2 / '['*3*']'",
                                                  "3 / '['*4*']'", 
                                                  "3 / '['*5*']'",
                                                  "4 / '['*5*']'"),], aes(x=mu,
                                                                                  y = TPR,
                                                                                  group = Methods,
                                                                                  linetype=Methods,
                                                                                  shape=Methods,
                                                                                  color=Methods)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(un ~ xcoef, labeller = label_parsed) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

pdf(file = "z u < n Results.pdf",
    width = 9.6, height = 9.6)
ggarrange(FDR_plot,TPR_plot,ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()