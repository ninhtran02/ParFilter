# Par-filter
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
obj <- parfilter(p = DMD.pvalues, error_targets = rep(0.05, 5), u = 3, groups = list(c(1,3),c(2,4)), 
                 group_options = c(2,1), selections = NULL, u_groups = c(2,1),
                 adaptive = TRUE, lambda = NULL,
                 auto = FALSE, omega = 0.5, w = c(0.5,0.5),
                 partition = NULL, method = "Stouffer")

# Print the results
print(obj)
```

## How to reproduce the simulation results for "Testing for Replicating Signals across Multiple Studies via Partioning and Filtering"
To reproduce the simulation results in an efficient manner, we assume the reader has access to an  account in a high performance computing (HPC) system running the *Slurm Workload Manager*. Follow the steps below:

1. Log onto your own HPC account.

2. Create a slurm file called *job_submission_bimodal.slurm* as follows:
```
#!/bin/bash
#SBATCH --job-name=REP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10000
#SBATCH --time=00-60:00:00

rho=$1
mu=$2
setup=$3

module --force purge
module load foss/2022a R/4.2.2

Rscript --vanilla Repmain.R $rho $mu $setup

```

Furthermore, create a slurm file called *batch_submission_bimodal.slurm* as follows:
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

3. In your HPC account, change your current directory to the ZDIRECT folder using the "cd" command:
```
cd    (your own working directory)/ZDIRECT
```
Make sure *main_bimodal2.R*, *job_submission_bimodal.slurm* and *batch_submission_bimodal.slurm* are inside the ZDIRECT folder.

4. To submit the jobs, run the two files *batch\_submission\_bimodal.slurm* with the command:
```
sbatch batch_submission_bimodal.slurm
```

5. The simulations will typically be finished in about a day, with the resulting data files saved to the subfolders within *ZDIRECT/simResults*, *ZDIRECT/simResultsrhopos05*, *ZDIRECT/simResultsrhopos08*, *ZDIRECT/simResultsrhoneg05* and *ZDIRECT/simResultsrhoneg08*.
