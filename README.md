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
obj <- ParFilter_FDR(p_mat = P, X_list = X_list, u = 4, q = 0.05, K = 4,
                             method = "Stouffer", adaptive = TRUE, cross_weights = FALSE,
                             lambdas = rep(0.50,K))

# Print the results
print(obj)
```
### Arguments
- `p_mat`: mxn matrix of p-values.
- `X_list`: list of length n, containing the covariates for each study.
- `u`: replicability threshold.
- `q`: FDR target.
- `K`: number of groups. ParFilter_FDR will automatically partition the n studies in two K groups of approximately equal sizes.
- `method`: Combining function for creating the local GBHPC p-values. Can be either: "Fisher", "Stouffer", or "Simes".
- `adaptive`:  logical indicating whether to use adaptive null proportion estimates or not.
- `cross_weights`: Set as TRUE if the p-values are dependent within studies, otherwise leave it as FALSE.
- `lambdas:` numeric of tuning parameters for the null proportion estimates.

### Values
`ParFilter_FDR` returns a numeric of features indices to be rejected.

## How to reproduce the simulation results for "Testing for Replicating Signals across Multiple Studies with Side Information"
To reproduce the simulation results in an efficient manner, we assume the reader has access to an account in a high performance computing (HPC) system running the *Slurm Workload Manager*. Follow the steps below:

1. Log onto your own HPC account.

2. Create a slurm file called *Repmain.slurm* as follows:
```
#!/bin/bash
#SBATCH --job-name=REP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5000
#SBATCH --time=00-24:00:00

xcoef=$1
mu=$2
u_n=$3

module --force purge
module load foss/2022a R/4.2.2

Rscript --vanilla Repmain.R $xcoef $mu $u_n
```
Feel free to edit the #SBATCH commands to suit your preferences. For example, you can insert `#SBATCH --mail-user=<your_email_address>` to send you a reminder email for when your simulation finishes.

Furthermore, create a slurm file called *batch_submission_Repmain.slurm* as follows:
```
#!/bin/bash

for rho in {1..3}
do

for mu in {1..6}
do

for setup in {1..9}
do

sbatch Repmain.slurm $rho $mu $setup

done

done

done

``` 

3. In your HPC account, change your current directory to the "Paper Simulations" folder using the "cd" command:
```
cd    (your own working directory)/Paper Simulations
```
Download *Repmain.R* in https://github.com/ninhtran02/ParFilter/tree/main/R and place it in the "Paper Simulations" folder. Furthermore, sake sure *main_bimodal2.R*, *Repmain.slurm* and *batch\_submission\_Repmain.slurm* are inside the folder as well.

4. To submit the jobs, run the file *batch\_submission\_Repmain.slurm* with the command:
```
sbatch batch_submission_Repmain.slurm
```

5. The simulations will typically be finished in about a day, with the resulting data files saved to the subfolders within *Paper Simulation/SavedData/Independence*, *Paper Simulation/SavedData/NegativeDependence*, and *Paper Simulation/SavedData/PositiveDependence*.
   
6. To run produce the plots, download and run *Plot.R*, *NegativeDependencPlot.R*, and *PositiveDependencePlot.R* from https://github.com/ninhtran02/ParFilter/tree/main/R.

## How to reproduce the real data results for "Testing for Replicating Signals across Multiple Studies via Partioning and Filtering"

Download the folder https://github.com/ninhtran02/ParFilter/tree/main/Real%20Data and run DMD_meta_analysis.R and AIRE_meta_analysis.R to reproduce the real data results.








