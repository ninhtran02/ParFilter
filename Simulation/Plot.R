setwd("~/Paper Simulations/")

library(ggplot2)
library(ggpubr)
library(latex2exp)

FDR_df <- data.frame(matrix(ncol = 4, nrow = 0))
names(FDR_df) <- c("FDR", "method", "mu", "(u,n)")
TPR_df <- data.frame(matrix(ncol = 4, nrow = 0))
names(TPR_df) <- c("TPR", "method", "mu", "(u,n)")

for(mu in c(2,2.3,2.6,2.9,3.2,3.5)){
  for(setup in list(c(2,2),c(3,4),c(4,4),
                    c(4,6),c(5,6),c(6,6),
                    c(6,8),c(7,8),c(8,8))){
    print(mu)
    u <- setup[1]
    n <- setup[2]
    file.name <- paste("SavedData/Independence/","mu",mu,"setup",u,n,".RD",sep = "")
    if(!file.exists(file.name)){
      FDR_dat <- cbind.data.frame(FDRrep,names(FDRrep), mu, paste("(",u,",",n,")", sep = "") )
      TPR_dat <- cbind.data.frame(TPRrep,names(TPRrep), mu, paste("(",u,",",n,")", sep = "") )
    } else {
      load(file = file.name)
      print(FDRrep)
      FDR_dat <- cbind.data.frame(FDRrep,names(FDRrep), mu, paste("(",u,",",n,")", sep = "") )
      TPR_dat <- cbind.data.frame(TPRrep,names(TPRrep), mu, paste("(",u,",",n,")", sep = "") )
    }
    FDRstudy
    names(FDR_dat) <- c("FDR", "method", "mu", "(u,n)")
    names(TPR_dat) <- c("TPR", "method", "mu", "(u,n)")
    
    FDR_df <- rbind.data.frame(FDR_df, FDR_dat)
    TPR_df <- rbind.data.frame(TPR_df, TPR_dat)
  }
}

rownames(FDR_df) <- NULL
rownames(TPR_df) <- NULL
methods <- c("PF", "MaxPF", "APF", "MaxAPF", "AutoAPF", "BH-Simes",
             "BH-Stouffer", "BH-Fisher", "AdaFilter")
FDR_df$method <- factor(x = FDR_df$method, levels = methods)
levels(FDR_df$method) <- c("ParFilter", "Max-ParFilter", "Adaptive-ParFilter",
                           "Adaptive-Max-ParFilter", "Adaptive-Auto-ParFilter",
                           "BH-Simes", "BH-Stouffer", "BH-Fisher", "AdaFilter")
#FDR_df$mu <- as.factor(FDR_dat$mu)
#FDR_df$`(u,n)` <- as.factor(FDR_dat$`(u,n)`)

TPR_df$method <- factor(x = TPR_df$method, levels = methods)
levels(TPR_df$method) <- c("ParFilter", "Max-ParFilter", "Adaptive-ParFilter",
                           "Adaptive-Max-ParFilter", "Adaptive-Auto-ParFilter",
                           "BH-Simes", "BH-Stouffer", "BH-Fisher", "AdaFilter")
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

pdf(file = "IndependentResultsNEW.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 15) #
ggarrange(FDR_plot_A, FDR_plot_B,FDR_plot_C,
          TPR_plot_A, TPR_plot_B,TPR_plot_C,
          ncol=1, nrow=6, common.legend = TRUE, legend="bottom")
dev.off()


##############
library(reshape2)

for(ss in list(c(2,2),c(3,4),c(4,4),c(4,6),c(5,6),c(6,6),c(6,8),c(7,8),c(8,8))){
  FDRstudy_df <- data.frame(matrix(ncol = 4, nrow = 0))
  TPRstudy_df <- data.frame(matrix(ncol = 4, nrow = 0))
  names(FDRstudy_df) <- c("Study", "method", "FDR", "mu")
  names(TPRstudy_df) <- c("Study", "method", "TPR", "mu")
  for(mu in c(2,2.3,2.6,2.9,3.2,3.5)){
    for(setup in list(c(ss[1],ss[2]))){
      #print(mu)
      u <- setup[1]
      n <- setup[2]
      file.name <- paste("SavedData/Independence/","mu",mu,"setup",u,n,".RD",sep = "")
      load(file = file.name)
      
      fdrstudy_dat <- cbind.data.frame(melt(FDRstudy),mu)
      tprstudy_dat <- cbind.data.frame(melt(TPRstudy),mu)
      names(fdrstudy_dat) <- c("Study", "method", "FDR", "mu")
      names(tprstudy_dat) <- c("Study", "method", "TPR", "mu")
      
      FDRstudy_df <- rbind.data.frame(FDRstudy_df, fdrstudy_dat)
      TPRstudy_df <- rbind.data.frame(TPRstudy_df, tprstudy_dat)
    }
  }

  methods <- c("ParFilter", "Max-ParFilter", "Adaptive-ParFilter",
               "Adaptive-Max-ParFilter","Adaptive-Auto-ParFilter",
               "BH-Simes", "BH-Stouffer", "BH-Fisher", "AdaFilter")
  PFmethods <- c("ParFilter", "Max-ParFilter", "Adaptive-ParFilter",
                 "Adaptive-Max-ParFilter","Adaptive-Auto-ParFilter")
  FDRstudy_df$Study <- as.factor(FDRstudy_df$Study)
  TPRstudy_df$Study <- as.factor(TPRstudy_df$Study)
  FDRstudy_df$method <- factor(FDRstudy_df$method)
  TPRstudy_df$method <- factor(TPRstudy_df$method)
  levels(FDRstudy_df$method) <- methods
  levels(TPRstudy_df$method) <- methods
FDR_study_plot <- ggplot(data=FDRstudy_df[FDRstudy_df$method %in% PFmethods,],
                                                     aes(x=mu,
                                                     y = FDR,
                                                     group = method,
                                                     linetype=method,
                                                     shape=method,
                                                     color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(PFmethods))) + 
  scale_color_manual(values = 10:(10+length(PFmethods)) ) +
  facet_grid(. ~ Study) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(FDR$_{Study}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.01, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

TPR_study_plot <- ggplot(data=TPRstudy_df[TPRstudy_df$method %in% PFmethods,],
                         aes(x=mu,
                             y = TPR,
                             group = method,
                             linetype=method,
                             shape=method,
                             color=method)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values=seq(0,length(PFmethods))) + 
  scale_color_manual(values = 10:(10+length(PFmethods)) ) +
  facet_grid(. ~ Study) +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  ylab(TeX(r'(TPR$_{Study}$)')) +
  xlab( TeX(r'($\mu$)') ) +
  theme_bw() +
  #geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

fna <- paste(ss[1],ss[2],"IndependentResultsNEW.pdf", sep = "")
pdf(file = fna,   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) #
print(ggarrange(FDR_study_plot, TPR_study_plot,
          ncol=1, nrow=2, common.legend = TRUE, legend="bottom"))
dev.off()

}

FDR_df_demo <- FDR_df[FDR_df$`(u,n)` %in% c("(2,2)","(4,4)","(6,6)","(8,8)") &
       FDR_df$method %in% c("BH-Stouffer","ParFilter") &
       FDR_df$mu == 3.5,]


TPR_df_demo <- TPR_df[TPR_df$`(u,n)` %in% c("(2,2)","(4,4)","(6,6)","(8,8)") &
                        TPR_df$method %in% c("BH-Stouffer","ParFilter") &
                        TPR_df$mu == 3.5,]


FDR_plot_demo <- ggplot(data=FDR_df_demo, aes(x= `(u,n)`,
                                                           y = FDR,
                                                           group = method,
                                                           linetype=method,
                                                           shape= method,
                                                           color=method)) +
  geom_line() +
  geom_point() +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  scale_shape_manual(values=c(1-1,7-1)) + 
  scale_color_manual(values = c(10,16) ) +
  scale_x_discrete(expand = c(0.01,0)) +
  ylab(TeX(r'(FDR$_{rep}$)')) +
  xlab( TeX(r'($(u,n)$)') ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype='dotted', size = 1) +
  theme(legend.title=element_blank())

TPR_plot_demo <- ggplot(data=TPR_df_demo, aes(x= `(u,n)`,
                                              y = TPR,
                                              group = method,
                                              linetype=method,
                                              shape= method,
                                              color=method)) +
  geom_line() +
  geom_point() +
  #scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10)) +
  scale_shape_manual(values=c(1-1,7-1)) + 
  scale_color_manual(values = c(10,16) ) +
  scale_x_discrete(expand = c(0.01,0)) +
  ylab(TeX(r'(TPR$_{rep}$)')) +
  xlab( TeX(r'($(u,n)$)') ) +
  theme_bw() +
  theme(legend.title=element_blank())

pdf(file = "DemoPlot.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 3) #
ggarrange(FDR_plot_demo, TPR_plot_demo,
          ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()
