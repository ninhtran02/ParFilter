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
Rejections <- ParFilter_FDR(p_mat = P, X_list = X_list, u = 4, q = 0.05, K = 4,
                             method = "Stouffer", adaptive = TRUE, cross_weights = FALSE,
                             lambdas = rep(0.50,K))

# Print the results
print(Rejections)
```
### Arguments
- `p_mat`: mxn matrix of p-values.
- `X_list`: list of length n, containing the covariates for each study.
- `u`: replicability threshold.
- `q`: FDR target.
- `K`: number of groups. ParFilter_FDR will automatically partition the n studies in two K groups of approximately equal sizes.
- `method`: combining function for creating the local GBHPC p-values. Can be either: "Fisher", "Stouffer", or "Simes".
- `adaptive`:  logical indicating whether to use adaptive null proportion estimates or not.
- `cross_weights`: set as TRUE if the p-values are dependent within studies, otherwise leave it as FALSE.
- `lambdas:` numeric of tuning parameters for the null proportion estimates.

### Values
`ParFilter_FDR` returns a numeric of features indices to be rejected.

## Reproducing the numerical results in "Testing for Replicating Signals across Multiple Studies with Side Information"

- Steps for reproducing the results in Section 4 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Section%204%20Simulations.md).
- Steps for reproducing the results in Section 5 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Section%205%20Case%20Study.md).
- Steps for reproducing the results in Appendix E can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20F.1%20Simulations.md)
- Steps for reproducing the results in Appendix F.1 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20F.1%20Simulations.md).
- Steps for reproducing the results in Appendix F.2 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20F.2%20Simulations.md).
- Steps for reproducing the results in Appendix F.3 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20F.3%20Simulations.md).
- Steps for reproducing the results in Appendix F.4 can be found [here](https://github.com/ninhtran02/ParFilter/blob/main/Appendix%20F.4%20Simulations.md).
