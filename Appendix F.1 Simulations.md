# Reproducing the simulation results in Appendix F.1
To reproduce the simulation results in an efficient manner, we assume the reader has access to an account in a high performance computing (HPC) system running the *Slurm Workload Manager*. Follow the steps below:

1. Log onto your own HPC account.

2. In your HPC account, change your to your working directory using the "cd" command:
```
cd    (your own working directory)/Paper Simulations/Covariate-Assisted
```
3. Create the following directories in your working directory:
```
mkdir -p SavedData/Independence/
```

3. Create a slurm file called `Repmain.slurm` as follows:
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
This slurm file will be used as a template for submitting a simulation job to the HPC under the parameter settings `xcoef`, `mu`, and `u_n`.

Feel free to edit the `#SBATCH` commands to suit your preferences. For example, you can insert `#SBATCH --mail-user=<your_email_address>` to send you a reminder email for when the simulation finishes. This may be useful since the simulations do take a while. 

3. Create a slurm file called `batch_submission_Repmain.slurm` as follows:
```
#!/bin/bash

for xcoef in {1..3}
do

for mu in {1..5}
do

for u_n in {1..8}
do

sbatch Repmain.slurm $xcoef $mu $u_n

done

done

done

```
This slurm file will be used to submit a simulation job to the HPC under every combination of parameter settings explored for `xcoef`, `mu`, and `u_n`. Yet again, feel free to add any `#SBATCH` commands to suit your preferences.


Download *Repmain.R* from [here](https://github.com/ninhtran02/ParFilter/tree/main/Simulation) and place it in the "Paper Simulations/Covariate-Assisted" folder. 

5. To submit the jobs, run the file *batch\_submission\_Repmain.slurm* with the command:
```
sbatch batch_submission_Repmain_z.slurm
```
Generally speaking, you can expect the simulations to finish within half a day or so.

6. The resulting data files saved to the subfolders within *Paper Simulation/Covariate-Assisted/SavedData/Independence*.
   
7. To produce the plots, download and run *Plot full.R* from [here](https://github.com/ninhtran02/ParFilter/tree/main/Simulation).





