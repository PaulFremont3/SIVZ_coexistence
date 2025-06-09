#!/bin/bash
#SBATCH -J SIVRZ_molar_model_phytotypes                           # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=10G                        # Memory per core
#SBATCH  --time=10:00:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_SIVRZ_molar_model_phytotypes-%j.out                        # Combined output and error messages file

python SIVRZ_molar_model_phytotypes.py $1 $2 $3 $4 $5 $6 $7 $8 $9

# $1 phytoplankton type;
# $2 % fail of infection;
# $3 growth rate ratio (resistant versus susceptible);
# $4 model type: SIVRZ or SIVRZ_hill
# $5 ocean type 0 upwelling or oligotrophic
# $6 dz2 22.4 or 0
# $7 dv2
# $8 param to explore: phir or epsr or lp_phir or lp_epsr
# $9 phi_ratio in the case param==lp (else can be anything), 0 = full resitance


