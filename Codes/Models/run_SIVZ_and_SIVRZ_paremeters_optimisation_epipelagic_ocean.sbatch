#!/bin/bash
#SBATCH -J SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean                           # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=6GB                        # Memory per core
#SBATCH  --time=00:45:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean-%j.out                        # Combined output and error messages file

python SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py $1 $2 $3 $4

# $1 phytoplankton type 0 1 2 or 3;
# $2 model type SIVZ, SIVZ_intra, SIVRZ, SIVRZ_intra
# $3 POM inhibition of phi
# $4 chunk
