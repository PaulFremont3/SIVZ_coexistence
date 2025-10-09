#!/bin/bash
#SBATCH -J SIVZ_molar_model_phytotypes                           # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=20G                        # Memory per core
#SBATCH  --time=10:00:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_SIVZ_molar_model_phytotypes-%j.out                        # Combined output and error messages file

 
python SIVZ_molar_model_phytotypes.py $1 $2 $3 $4 $5 $6 $7 $8 $9

# $1 phytoplankton type 0 1 2 or 3;
# $2 loss of infectivity;
# $3 ocean type 0 oligotrophic or upwelling
# $4 quadratic mortalitiy of zoop (22.4 or 0)
# $5 quadratic mortalitiy of virus 
# $6 strength of intraguild predation 0-1
# $7 predation type: 0 =DARWIN, encounter_model=Talmy 2019
# $8 0 or 1; if 0=> smallest cell if 1 => biggest cell
# $9 epsilon



