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

## How to reproduce the simulation results for "TESTING FOR REPLICATING SIGNALS ACROSS MULTIPLE STUDIES VIA PARTITIONING AND FILTERING"
To reproduce the simulation results in an efficient manner, we assume the reader has access to an  account in a high performance computing (HPC) system running the *Slurm Workload Manager*. Follow the steps below:

1. Log onto your own HPC account.

2. Create a slurm file called *job_submission_bimodal.slurm* as follows:
```
#!/bin/bash
#SBATCH --job-name=dFDR_bimodal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3500
#SBATCH --time=00-80:00:00

nullprop=$1
mu=$2
symm=$3
rho=$4


module --force purge
module load foss/2022a R/4.2.2


Rscript --vanilla main_bimodal2.R $nullprop $mu $symm $rho

```

Furthermore, create a slurm file called *batch_submission_bimodal.slurm* as follows:
```
#!/bin/bash


for nullprop in 0 0.2 0.5 0.8
do

  for mu in 0.5 1 1.5 2 2.5 3 
  do
  
    for symm in 0.5 0.75 1
    do
    	for rho in -0.8	-0.5 0 0.50 0.80
    	do
  
              sbatch /data/gpfs/projects/punim1426/ZDIRECT/job_submission_bimodal.slurm $nullprop $mu $symm $rho
	
  	done
    
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
