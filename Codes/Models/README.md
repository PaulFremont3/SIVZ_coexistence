In this directory: all python codes to run simulations of the SVZ, SIVZ, SVRZ and SIVRZ models and the parameter otpimization procedure for 4 types of phytoplankton: *Prochlorococcus*, *Synechococcus*, a picoeukaryote, a small diatom.  

All simulations were run on the Zaratan HPC cluster of the University of Maryland (https://hpcc.umd.edu/hpcc/zaratan.html)

# 1. DESCRIPTION OF FILES  

## 1.1 Python files:
in each case pdf outputs are specified and *\*suffix\** refers to a suffix to identify parameters input of the simulation (see inside .py file and .sbatch files for inputs needed for each files): 
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

## 1.2 .sbatch files:
(each file contains a description of the parameters taken as inputs): to run simulations on HPC cluster
  - **run_SVZ_molar_model_phytotypes.sbatch**: runs SVZ_molar_model_phytotypes.py
  - **run_SVRZ_molar_model_phytotypes.sbatch**: runs SVRZ_molar_model_phytotypes.py
  - **run_SIVZ_molar_model_phytotypes.sbatch**: runs SIVZ_molar_model_phytotypes.py
  - **run_SIVRZ_molar_model_phytotypes.sbatch**: runs SIVRZ_molar_model_phytotypes.py
  - **run_SIVZ_MCT.sbatch** runs SIVZ_MCT.py
  - **run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch**: runs SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py
  - **run_analysis_optimization.sbatch**: runs analysis_optimization.py

## 1.3 .sh files: 
  - **run_coexistence_simulations.sh** : runs all simulations necessary to generate figure 2, 4, 6b-j, S5, S7, s8, s9, s10, s11, s12, s13, s17 and data for: figure 5, 3d,e, S6, s14, s16
  - **run_simulations_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sh**: used to run in parallel  the optimization of parameters, it divides the grid search in 50 equal chunks of parameters combinations e.g. if 10^6 parameter combination are tested, one job will test 10^6/50= 20000 combinations 
  - **concatenate_results_optimisation_params.sh**: used to concatenate rsult files from the optimisation

# 2. HOW TO RUN THE SIMULATIONS

In a HPC cluster:  

## 2.1. Create in this directory a subfolder called `model_data` to store results:  
&nbsp;&nbsp;&nbsp;&nbsp;`mkdir model_data`

## 2.2. Run simulations of dynamics across the parameter space:  
&nbsp;&nbsp;&nbsp;&nbsp;`./run_coexistence_simulations.sh`  
&nbsp;&nbsp;&nbsp;&nbsp;This step will run all simulations necessary to generate figures of the coexistence diagram and some data that needs to be rescaled to generate some figures.

&nbsp;&nbsp;&nbsp;&nbsp;**Run times:**  
&nbsp;&nbsp;&nbsp;&nbsp;- SVZ and SVRZ simulations across the adsorption rate parameter space: ~3 minutes  
&nbsp;&nbsp;&nbsp;&nbsp;- SIVZ and SIVRZ simulations across adsorption rate and latent period: ~2h30 to 3 hours  
&nbsp;&nbsp;&nbsp;&nbsp;- SVRZ and SIVRZ across adsorption rate and resistance strength: ~3h30  

&nbsp;&nbsp;&nbsp;&nbsp;Before running all together, consider checking each job (SVZ, SIVZ, SVRZ, and SIVRZ) individually to ensure it works.

## 2.4. Generate scaled figures:  
&nbsp;&nbsp;&nbsp;&nbsp;`python make_figures_scaled.py main_models`  
&nbsp;&nbsp;&nbsp;&nbsp;This will generate figures that need to be scaled (plots saved in the `model_data/` folder).

## 2.5. Create a subfolder to store parameter optimization results:  
&nbsp;&nbsp;&nbsp;&nbsp;`mkdir results_optimization_params/`

## 2.6. Grid search to optimize parameters to target concentrations:  
&nbsp;&nbsp;&nbsp;&nbsp;`./run_simulations_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sh 0`  
&nbsp;&nbsp;&nbsp;&nbsp;This will run optimization of parameters to minimize distance to target concentrations for the four phytoplankton types considered in the study:  
&nbsp;&nbsp;&nbsp;&nbsp;a small diatom, a picoeukaryote, a *Synechococcus*, and a *Prochlorococcus*.  
&nbsp;&nbsp;&nbsp;&nbsp;Results are stored in the `results_optimization_params/` folder.

## 2.7. Concatenate files of optimization results:  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Synechococcus 0`  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Prochlorochoccus 0`  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Eukaryote 0`  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Diatom 0`

## 2.8. Run the analysis of the optimization:  
&nbsp;&nbsp;&nbsp;&nbsp;`sbatch run_analysis_optimization.sbatch 0`  
&nbsp;&nbsp;&nbsp;&nbsp; This will extract the best and 200 best parameter combinations for each phytoplankton type

# 3. PDF AND TXT FILES WITH THE FIGURES AND DATA FOR TABLES IN THE PAPER 

(after running all simulations)  

# Figures
- Figure 2
  - 2a: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_no-dz2_mesotrophic.pdf
  - 2b: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_mesotrophic.pdf
  - 2c: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_no-dz2_m2-1800.0_mesotrophic.pdf
  - 2d: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 2e: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_no-dz2_mesotrophic.pdf
  - 2f: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_mesotrophic.pdf
  - 2g: page 1 of SIVZ_m2_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_m2-1800.0_no-dz2_mesotrophic.pdf
  - 2h: page 1 of SIVZ_m2_model_phi_latent_period_Prochlorochoccus_BS20.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 2i: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-10.0.pdf
  - 2j: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_mesotrophic_phir-10.0.pdf
  - 2k: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_phir-10.0.pdf
  - 2l: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_mesotrophic_phir-10.0.pdf
  - 2m: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_no-dz2_lp_phi-ratio-10.0_mesotrophic.pdf
  - 2n: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_lp_phi-ratio-10.0_mesotrophic.pdf
  - 2o: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_no-dz2_m2-1800_lp_phi-ratio-10.0_mesotrophic.pdf
  - 2p: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-10.0_mesotrophic.pdf
  - 2q: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-0.0.pdf
  - 2r: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_mesotrophic_phir-0.0.pdf
  - 2s: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_phir-0.0.pdf
  - 2t: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS20.0_LOI0_GR-R0.8_mesotrophic_phir-0.0.pdf
  - 2u: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_no-dz2_lp_phi-ratio-0.0_mesotrophic.pdf
  - 2v: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_lp_phi-ratio-0.0_mesotrophic.pdf
  - 2w: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_no-dz2_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf
  - 2x: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS20.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf

