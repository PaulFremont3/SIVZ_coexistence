#!/usr/bin/env bash
#SBATCH -J analysis_optimization                           # Job name
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=30G                        # Memory per core
#SBATCH  --time=00:10:00                        # Duration of the job (Ex: 10 mins)
#SBATCH -o run_analysis_optimization-%j.out                        # Combined output and error messages file

echo "Job started at $(date)"

python analysis_optimization.py $1

echo "Job finished at $(date)"
