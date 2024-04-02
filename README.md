# Par-filter
 Implementation of the PC partitioning and filtering algorithm

 ## EXAMPLE
 ```
install.packages("devtools")
devtools::install_github("ninhtran02/Par-filter")
 obj <- Repfilter(p = p.mat, error_targets = rep(0.05, 5), u = 3, groups = list(c(1,3),c(2,4)), 
                  group_options = c(2,1), selections = NULL, u_groups = c(2,1),
                  adaptive = TRUE, lambda = NULL,
                  auto = FALSE, omega = 0.75, w = c(0.5,0.5),
                  partition = NULL, method = "Simes")
print(obj)
```
