# ParFilter
 Implementation of the partial conjunction partitioning and filtering algorithm.

 ## Installation
```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("ninhtran02/Parfilter")
```

 ## Usage example
 ```
library(ParFilter)
m <- 5000
n <- 4
P <- matrix(c(rnorm(19000,0),rnorm(1000,3)), nrow = m, ncol = n, byrow = TRUE)
P <- 1 - pnorm(P)
X_list <- lapply(1:4, function(x) rnorm(5000)) # Non-informative covariates
Rejections <- ParFilter(P = P, u = 4, X_list = X_list, K = 4, q = 0.05, lambda_vec = rep(0.50,4))

# Print the results
print(Rejections)
```
### Arguments
- `P`: mxn matrix of p-values.
- `u`: replicability threshold.
- `X_list`: list of length n, containing the covariates for each study.
- `K`: number of groups. ParFilter_FDR will automatically partition the n studies in two K groups of approximately equal sizes.
- `direction`: If u < n and covariates are univariate, set `direction = negative` if the smaller the covariate, the more likely its corresponding base null is false, and set `direction = positive` if the larger the covariate, the more likely its corresponding base null is false.
- `q`: FDR target.
- `inflate`:  logical indicating whether to inflate the null proportion estimate to account for dependence among base p-values.
- `lambda_vec:` numeric of tuning parameters for the null proportion estimates.

### Values
`ParFilter_FDR` returns a numeric of features indices to be rejected.

## Reproducing the numerical results in "Testing for Replicating Signals across Multiple Studies with Side Information"

- Steps for reproducing the results in Section 4 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Section%204%20Simulations.md).
- Steps for reproducing the results in Section 5 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Section%205%20Case%20Study.md).
- Steps for reproducing the results in Appendix D can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20D%20simulations.md)
- Steps for reproducing the results in Appendix E.1 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20E.1%20Simulations.md).
- Steps for reproducing the results in Appendix E.2 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20E.2%20Simulations.md).
- Steps for reproducing the results in Appendix E.3 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20E.3%20Simulations.md).
- Steps for reproducing the results in Appendix E.4 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20E.4%20Simulations.md).
- Steps for reproducing the results in Appendix E.5 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20E.4%20Simulations.md).
