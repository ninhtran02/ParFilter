Test

```
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

demo_plot <- function(n, n_sims, m, q){
  m_vec <- c(1,seq(10,m,10))
  u_vec <- 1:n
  
  dat_mat <- matrix(data = 0, nrow = length(m_vec), ncol = length(u_vec))
  
  for(iter in 1:n_sims){
    print(iter)
    for(m in m_vec){
      for(u in u_vec){
        pc_p_tracker <- c()
        for(i in 1:m){
          p_vec <- rbeta(n = n, shape1 = 0.26, shape2 = 7)
          # Change combining function here
          # ParFilter::fisher_fun(p_vec = p_vec, u = u) for Fisher's methhod
          # ParFilter::stouffer_fun(p_vec = p_vec, u = u) for Stouffer's method
          pc_p <- ParFilter::fisher_fun(p_vec = p_vec, u = u)
          pc_p_tracker <- c(pc_p_tracker,pc_p)
        }
        rej <- sum(p.adjust(p = pc_p_tracker, method = "BH") <= q)/m
        dat_mat[which(m_vec == m), which(u_vec == u)] <- dat_mat[which(m_vec == m), which(u_vec == u)] + rej/n_sims
      }
    }
  }
  
  dat_df <- as.data.frame(dat_mat)
  colnames(dat_df) <- 1:u
  dat_df <- reshape2::melt(dat_df)
  dat_df$m <- rep(m_vec, length(u_vec))
  dat_df$variable <- as.factor(dat_df$variable)
  dat_df$n <- rep(paste("n =",n),nrow(dat_df))
  colnames(dat_df) <- c("u","Power","m","n")
  
  return(dat_df)
}
 
n_sims <- 500

m <- 1000
q <- 0.05

N <- 8
overall_dat <- data.frame()
plot_list <- rep(list(),N)
for(n in 2:N){
  overall_dat <- rbind(overall_dat, demo_plot(n = n, n_sims = n_sims, m = m, q = 0.05))
}

PCdemo <- ggplot(overall_dat, aes(x = m, y = Power)) + 
  geom_line(aes(color = u, linetype = u), size = 0.75) +
  ylab(TeX(r'(Power)')) +
  xlab(TeX(r'($m$)')) +
  facet_wrap(~n, ncol = 1, strip.position="right") +
  theme_bw() +
  theme(legend.position = "bottom") +             # position at bottom
  scale_x_continuous(
    breaks = c(1, 250, 500, 750, 1000),
    limits = c(1, 1000),
    expand = c(0, 0)
  ) + 
  guides(colour = guide_legend(nrow = 1))

pdf(file = "PCdemos_more_n",
    width = 9.6, height = 9.6)
PCdemo
dev.off()
```
