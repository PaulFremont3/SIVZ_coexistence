#!/bin/bash
#SBATCH -J SIVZ_MCT                          # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=5G                        # Memory per core
#SBATCH  --time=24:00:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_SIVZ_MCT-%j.out                        # Combined output and error messages file
 
python SIVZ_MCT.py $1 $2 $3 $4 $5 $6

# $1 phytoplankton type 0 1 2 or 3;
# $2 loss of infectivity;
# $3 ocean type 0 oligotrophic, mesotrophic or upwelling
# $4 dz2 0 or 22.4
# $5 strength of intraguild predation 0-1
# $6 dv2



