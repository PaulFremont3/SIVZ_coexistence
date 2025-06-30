#!/bin/bash
#SBATCH -J ml_model                           # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=5G                        # Memory per core
#SBATCH  --time=24:00:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_ml_model-%j.out                        # Combined output and error messages file

module load r/4.2.1-bio 
Rscript best_model_ml.R $1 $2




