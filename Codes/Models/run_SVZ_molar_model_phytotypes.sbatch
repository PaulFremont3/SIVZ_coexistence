#!/bin/bash
#SBATCH -J SVZ_molar_model_phytotypes                           # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=5G                        # Memory per core
#SBATCH  --time=06:00:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_SVZ_molar_model_phytotypes-%j.out                        # Combined output and error messages file

echo "Job started at $(date)"

python SVZ_molar_model_phytotypes.py $1 $2 $3 $4 $5

echo "Job finished at $(date)"

# $1 phytoplankton type 0 1 2 or 3;
# $2 loss of infectivity;
# $3 ocean type 0 oligotrophic or upwelling
# $4 quadratic mortalitiy of zoop (22.4 or 0)
# $5 quadratic mortality of virus



