In this directory: all python codes to run simulations of the SVZ, SIVZ, SVRZ and SIVRZ models

**DESCRIPTION OF FILES**
- Python files, in each case pdf outputs are specified and *\*suffix\** refers to a suffix to identify parameters input of the simulation (see inside .py file and .sbatch files for inputs needed for each files): 
  - **generic_functions.py**: files containing functions shared by different scripts
  - **SVZ_molar_model_phytotypes.py** : code to run SVZ simulations and generate coexistence diagrams
      - outputs:
        - SVZ_model_phi_latent_period_time_series_*\*suffix\**.pdf: time series for different adsorption rates
        - SVZ_model_phi_latent_period_*\*suffix\**.pdf: coexistence diagram and other features in the parameter space
  - **SVZ_functions.py** : functions specific to the SVZ model
  - **SVRZ_molar_model_phytotypes.py** : code to run SVRZ simulations and generate coexistence diagrams
    - outputs:
      - if dz2==0, SVRZ_model_phi_latent_period_time_series_coex_*\*suffix\**.pdf: if the quadratic mortality term of the zooplankton is 0: time series in the case case of coexistence
      - SVRZ_model_phi_latent_period_time_series_*\*suffix\**.pdf: time series for different adsorption rates with phi/phir=10
      - SVRZ_model_phi_latent_period_time_series_full-res_*\*suffix\**.pdf: time series for different adsorption rates with phir=10 (full resistance)
      - SVRZ_model_phi_versus_phir_*\*suffix\**.pdf: coexistence diagram and other features in the parameter space
  - **SVRZ_functions.py** : functions specific to the SVRZ model
  - **SIVZ_molar_model_phytotypes.py** : code to run SIVZ simulations and generate coexistence diagrams
    - outputs:
      - if dz2==dv2==0, SIVZ_model_phi_latent_period_time_series_coex_*\*suffix\**.pdf: if the quadratic mortality terms are 0, time series in the case case of coexistence
      - SIVZ_model_phi_latent_period_time_series_*\*suffix\**.pdf: time series for different adsorption rates, fixed latent period (life history trait model)
      - SIVZ_model_phi_latent_period_*\*suffix\**.pdf: coexistence diagram and other features in the parameter space
  - **SIVZ_functions.py** : functions specific to the SIVZ model
  - **SIVRZ_molar_model_phytotypes.py** : code to run SIVRZ simulations and generate coexistence diagrams
    - outputs:
      - SIVRZ_model_phi_latent_period_time_series_coex_*\*suffix\**.pdf: if the quadratic mortality terms are 0, time series in the case case of coexistence
      - if param=='phir' or param=='epsr' (param is an input parameter), SIVRZ_model_phi_latent_period_time_series_*\*suffix\**.pdf, in case the resistance strength parameter space is explored, time series for different adsorption rates with phi/phir=10
      - if param=='phir' or param=='epsr', SIVRZ_model_phi_latent_period_time_series_full-res_*\*suffix\**.pdf', in case the resistance strength parameter space is explored, time series for different adsorption rates with phir=0 (full resistance)
      - Main pdf:
        - if param=='phir', SIVRZ_model_phi_versus_phir_*\*suffix\**.pdf, in the case extracellular resistance space is explored, coexistence diagram and other features in the parameter space
        - if param=='epsr', SIVRZ_model_eps_versus_epsr_*\*suffix\**.pdf, in the case extracellular resistance space is explored, coexistence diagram and other features in the parameter space
        - if param=='lp_phir' or param=='lp_epsr', in the case latent period space is explored (fixed intra or extra cellular resistance), coexistence diagram and other features in the parameter space
  Note: the SVZ_molar_model_phytotypes.py, SVRZ_molar_model_phytotypes.py SIVZ_molar_model_phytotypes.py, and SIVRZ_molar_model_phytotypes.py also save as .txt files some of the matrices they generate
  - **SIVRZ_functions.py** : functions specific to the SIVRZ model
  - **SIVZ_MCT.py**: code to run the SIVZ Modern coexistence theory analysis
    - outputs:
      - SIVZ_model_phi_latent_period_*\*suffix\**_MCT.pdf, effects from MCT Theory across the parameter space
  - **make_figures_scaled.py**: code to create figures with same scales (from different models): => figure 5, s6, s14, and s16
    - outputs:
      - Scaled_PVZ_main_models.pdf
      - Scaled_PVZ_FFT_and_Perdiod.pdf
      - Scaled_SIVZ_coexistence.pdf
      - Scaled_SIVZ_distances.pdf
  - **SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py**: file to run the optimization of parameterization of the SIVZ and SIVRZ models with respect to target concentrations (grid search)
  - **analysis_optimization.py**: analyze the results of the optimization =>  generates data for table 1 and 2

- .sbatch files (each file contains a description of the parameters taken as inputs): to run simulations on HPC cluster
  - **run_SVZ_molar_model_phytotypes.sbatch**: runs SVZ_molar_model_phytotypes.py
  - **run_SVRZ_molar_model_phytotypes.sbatch**: runs SVRZ_molar_model_phytotypes.py
  - **run_SIVZ_molar_model_phytotypes.sbatch**: runs SIVZ_molar_model_phytotypes.py
  - **run_SIVRZ_molar_model_phytotypes.sbatch**: runs SIVRZ_molar_model_phytotypes.py
  - **run_SIVZ_MCT.sbatch** runs SIVZ_MCT.py
  - **run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch**: runs SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py
  - **run_analysis_optimization.sbatch**: runs analysis_optimization.py

- .sh files: 
  - **run_coexistence_simulations.sh** : runs all simulations necessary to generate figure 2, 4, 6b-j, S5, S7, s8, s9, s10, s11, s12, s13, s17 and data for: figure 5, 3d,e, S6, s14, s16
  - **run_simulations_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sh**: used to run in parallel  the optimization of parameters, it divides the grid search in 50 equal chunks of parameters combinations e.g. if 10^6 parameter combination ar tested, one job will test 10^6/50= 20000 combinations 
  - **concatenate_results_optimisation_params.sh**: used to concatenate rsult files from the optimisation

**HOW TO RUN THE SIMULATIONS**
1. Create in this directorry a subfolfder called model_data to store results: `mkdir model_data`
2. run simulations of dynamics across the parameter space: command `./run_coexistence_simulations.sh`
   This step will run all simulation necessary to generate figures of coexistence diagram an some data that needs to be rescaled to generate some figures  
   Run times  
   SVZ and SVRZ simulations across the adsorption rate parameter space ~ 3 minutes  
   SIVZ and SIVRZ simulations across the adsorption rate and latent period parameter space ~ 2h30 to 3 hours  
   SVRZ and SIVRZ across the adsorption rate and resistance strenght parameter space ~ 3 hour 30  
   Before running all together, eventually check individually each type of job (SVZ, SIVZ, SVRZ and SIVRZ) to see whether it works  
4. Generate scaled figures: `python make_figures_scaled.py main_models`
   This will generate figures that need to be scaled (plots saved data from the model_data/ folder)  
5. Create in this directory a subfolder called results_optimization_params to store parameter optimization results: `mkdir results_optimization_params/`  
6. Grid search to optimize parameters to target concentrations: `./run_simulations_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sh 0`  
  This will run the optimization of parameters to minimize distance to target concentration for the four phytoplankton types considered in the study: a small diatom, a picoeukaryote, a synechococcus and a prochlorococcus  
7.  

**WHERE ARE THE FIGURES FROM THE PAPER**
