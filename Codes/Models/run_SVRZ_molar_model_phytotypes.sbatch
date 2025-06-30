#!/bin/bash
#SBATCH -J SVRZ_molar_model_phytotypes                           # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=12G                        # Memory per core
#SBATCH  --time=10:00:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_SVRZ_molar_model_phytotypes-%j.out                        # Combined output and error messages file

echo "Job started at $(date)"

python SVRZ_molar_model_phytotypes.py $1 $2 $3 $4 $5 $6 $7

echo "Job finished at $(date)"

# $1 phytoplankton type 0 1 2 or 3;
# $2 ocean type 0 oligotrophic or upwelling;
# $3 resistant growth rate penalty
# $4 quadratic mortalitiy of zoop (22.4 or 0)
# $5 quadratic mortalitiy of virus
# $6 resistance type: phir or epsr
# $7 specify a specific value of phi/phir or eps/epsr => plot 1d matrix, If 'no' => explore the phi/phir or eps/epsr space



