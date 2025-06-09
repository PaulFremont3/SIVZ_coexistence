In this directory: all python codes to run simulations of the SVZ, SIVZ, SVRZ and SIVRZ models

- Python files: 
  - generic_functions.py: files containing functions shared by different scripts
  - SVZ_molar_model_phytotypes.py : code to run SVZ simulations and generate coexistence diagrams
  - SVZ_functions.py : functions specific to the SVZ model
  - SVRZ_molar_model_phytotypes.py : code to run SVRZ simulations and generate coexistence diagrams
  - SVRZ_functions.py : functions specific to the SVRZ model
  - SIVZ_molar_model_phytotypes.py : code to run SIVZ simulations and generate coexistence diagrams
  - SIVZ_functions.py : functions specific to the SIVZ model
  - SIVRZ_molar_model_phytotypes.py : code to run SIVRZ simulations and generate coexistence diagrams
  - SIVRZ_functions.py : functions specific to the SIVRZ model
  - SIVZ_MCT.py: code to run the SIVZ Modern coexistence theory analysis
  - make_figures_scaled.py: code to create figures with comparable scales (from different models): => figure 5, s6, and s14
  - SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py: file to run the optimization of parameterization of the SIVZ and SIVRZ models with respect to target concentrations (grid search)
  - analysis_optimization.py: analyze the results of the optimization =>  generates data for table 1 and 2

- .sbatch files (each file contains a description of the parameters taken as inputs): to run simulations on HPC cluster
  - run_SVZ_molar_model_phytotypes.sbatch: runs SVZ_molar_model_phytotypes.py
  - run_SVRZ_molar_model_phytotypes.sbatch: runs SVRZ_molar_model_phytotypes.py
  - run_SIVZ_molar_model_phytotypes.sbatch: runs SIVZ_molar_model_phytotypes.py
  - run_SIVRZ_molar_model_phytotypes.sbatch: runs SIVRZ_molar_model_phytotypes.py
  - run_SIVZ_MCT.sbatch runs SIVZ_MCT.py
  - run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch: runs SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py
  - run_analysis_optimization.sbatch: runs analysis_optimization.py

- .sh files: 
  - run_coexistence_simulations.sh : runs all simulations necessary to generate figure 2, 3d,e, 4, 6b-j, S5, S7, s8, s9, s10, s11, s12, s13, s16, s17 and data for: figure 5, S6, s14
  - run_simulations_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sh: used to run in parallel  the optimization of parameters, it divides the grid search in 50 equal chunks of parameters combinations e.g. if 10^6 parameter combination ar tested, one job will test 10^6/50= 20000 combinations 
  - concatenate_results_optimisation_params.sh: used to concatenate rsult files from the optimisation
