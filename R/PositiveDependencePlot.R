setwd("~/Paper Simulations/")

library(ggplot2)
library(ggpubr)
library(latex2exp)

FDR_df <- data.frame(matrix(ncol = 4, nrow = 0))
names(FDR_df) <- c("FDR", "method", "mu", "(u,n)")
TPR_df <- data.frame(matrix(ncol = 4, nrow = 0))
names(TPR_df) <- c("TPR", "method", "mu", "(u,n)")

methods <- c("PF", "MaxPF", "APF", "MaxAPF", "AutoAPF",
             "BY-Simes", "BY-Stouffer", "BY-Fisher", "AdaFilter", "ad-PF", "prds-PF")

for(mu in c(2,2.3,2.6,2.9,3.2,3.5)){
  for(setup in list(c(2,2),c(3,4),c(4,4),
                    c(4,6),c(5,6),c(6,6),
                    c(6,8),c(7,8),c(8,8))){
    print(mu)
    u <- setup[1]
    n <- setup[2]
    file.name <- paste("SavedData/PositiveDependence/","mu",mu,"setup",u,n,".RD",sep = "")
    if(!file.exists(file.name)){
      FDR_dat <- cbind.data.frame(FDRrep,names(FDRrep), mu, paste("(",u,",",n,")", sep = "") )
      TPR_dat <- cbind.data.frame(TPRrep,names(TPRrep), mu, paste("(",u,",",n,")", sep = "") )
    } else {
      load(file = file.name)
      #print(FDRrep)
      FDR_dat <- cbind.data.frame(FDRrep,names(FDRrep), mu, paste("(",u,",",n,")", sep = "") )
      TPR_dat <- cbind.data.frame(TPRrep,names(TPRrep), mu, paste("(",u,",",n,")", sep = "") )
    }
    names(FDR_dat) <- c("FDR", "method", "mu", "(u,n)")
    names(TPR_dat) <- c("TPR", "method", "mu", "(u,n)")
    FDR_dat <- FDR_dat[FDR_dat$method %in% methods,]
    TPR_dat <- TPR_dat[TPR_dat$method %in% methods,]
    
    FDR_df <- rbind.data.frame(FDR_df, FDR_dat)
    TPR_df <- rbind.data.frame(TPR_df, TPR_dat)
  }
}

rownames(FDR_df) <- NULL
rownames(TPR_df) <- NULL
FDR_df$method <- factor(x = FDR_df$method, levels = methods)
#FDR_df$mu <- as.factor(FDR_dat$mu)
#FDR_df$`(u,n)` <- as.factor(FDR_dat$`(u,n)`)

TPR_df$method <- factor(x = TPR_df$method, levels = methods)
#TPR_df$mu <- as.factor(TPR_dat$mu)
#TPR_df$`(u,n)` <- as.factor(TPR_dat$`(u,n)`)

rownames(FDR_df) <- NULL
rownames(TPR_df) <- NULL
#methods <- c("PF", "MaxPF", "APF", "MaxAPF", "AutoAPF", "BH-Simes",
#             "BH-Stouffer", "BH-Fisher", "AdaFilter")
FDR_df$method <- factor(x = FDR_df$method)
levels(FDR_df$method) <- c("ParFilter", "Max-ParFilter", "Adaptive-ParFilter",
                           "Adaptive-Max-ParFilter", "Adaptive-Auto-ParFilter",
                           "BY-Simes", "BY-Stouffer", "BY-Fisher", "AdaFilter",
                           "ArbDep-ParFilter", "PRDS-ParFilter")
#FDR_df$mu <- as.factor(FDR_dat$mu)
#FDR_df$`(u,n)` <- as.factor(FDR_dat$`(u,n)`)

TPR_df$method <- factor(x = TPR_df$method, levels = methods)
#TPR_df$mu <- as.factor(TPR_dat$mu)
#TPR_df$`(u,n)` <- as.factor(TPR_dat$`(u,n)`)

FDR_plot_A <- ggplot(data=FDR_df[FDR_df$`(u,n)` %in% c("(2,2)","(3,4)","(4,4)") ,], aes(x=mu,
                                                                                        y = FDR,
                                                                                        group = method,
                                                                                        linetype=method,
                                                                                        shape=method,
                                                                                        color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(. ~ `(u,n)`) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

FDR_plot_B <- ggplot(data=FDR_df[FDR_df$`(u,n)` %in% c("(4,6)","(5,6)","(6,6)"),], aes(x=mu,
                                                                                       y = FDR,
                                                                                       group = method,
                                                                                       linetype=method,
                                                                                       shape=method,
                                                                                       color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(. ~ `(u,n)`) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

FDR_plot_C <- ggplot(data=FDR_df[FDR_df$`(u,n)` %in% c("(6,8)","(7,8)","(8,8)"),], aes(x=mu,
                                                                                       y = FDR,
                                                                                       group = method,
                                                                                       linetype=method,
                                                                                       shape=method,
                                                                                       color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(. ~ `(u,n)`) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())



ggarrange(FDR_plot_A, FDR_plot_B, FDR_plot_C, ncol=1, nrow=3, common.legend = TRUE, legend="bottom")

TPR_plot_A <- ggplot(data=TPR_df[TPR_df$`(u,n)` %in% c("(2,2)","(3,4)","(4,4)"), ], aes(x=mu,
                                                                                        y = TPR,
                                                                                        group = method,
                                                                                        linetype=method,
                                                                                        shape=method,
                                                                                        color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(. ~ `(u,n)`) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab( TeX(r'(TPR$_{rep}$)') ) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

TPR_plot_B <- ggplot(data=TPR_df[TPR_df$`(u,n)` %in% c("(4,6)","(5,6)","(6,6)"), ], aes(x=mu,
                                                                                        y = TPR,
                                                                                        group = method,
                                                                                        linetype=method,
                                                                                        shape=method,
                                                                                        color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(. ~ `(u,n)`) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab( TeX(r'(TPR$_{rep}$)') ) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

TPR_plot_C <- ggplot(data=TPR_df[TPR_df$`(u,n)` %in% c("(6,8)","(7,8)","(8,8)"), ], aes(x=mu,
                                                                                        y = TPR,
                                                                                        group = method,
                                                                                        linetype=method,
                                                                                        shape=method,
                                                                                        color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(methods))) + 
  scale_color_manual(values = 10:(10+length(methods)) ) +
  facet_grid(. ~ `(u,n)`) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab( TeX(r'(TPR$_{rep}$)') ) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

pdf(file = "PositiveDependenceResultsNEW.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 15) #
ggarrange(FDR_plot_A, FDR_plot_B,FDR_plot_C,
          TPR_plot_A, TPR_plot_B,TPR_plot_C,
          ncol=1, nrow=6, common.legend = TRUE, legend="bottom")
dev.off()
