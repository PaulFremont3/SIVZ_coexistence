In this directory: all python and R codes to:
- run simulations of the SVZ, SIVZ, SVRZ and SIVRZ models,
- example data in the folder `model_data/` from Figure 4, Figure 7a-o, Figure S7, Figure S9, Figure S10, Figure S12, Figure S14, and Figure S16 are stored in the folder `model_data` so these figures can be directly reproduced by running:  
   - `python make_figures_scaled.py main_models` for Figure 7a-o, Figure S7, Figure S14 and Figure S16: generating the files `Scaled_PVZ_FFT_and_Period.pdf`,  `Scaled_PVZ_main_models.pdf`,  `Scaled_SIVZ_coexistence.pdf`, and  `Scaled_SIVZ_distances.pdf` in <10 seconds)
   -  `python make_state_diagrams_figures.py` for Figure 4, Figure S9, Figure S12 and Figure S14: generating the files `State_diagrams_main_models.pdf`, `State_diagrams_resistance_models.pdf`, `State_diagrams_phytoplankton_types_models.pdf`, `State_diagrams_phytoplankton_types_resistance_models.pdf`, `State_diagrams_phytoplankton_types_th.pdf`, and `State_diagrams_phytoplankton_types_resistance_th.pdf`)  
- to reproduce other figures, simulations will have to be run
- the parameter otpimization procedure for 4 types of phytoplankton: *Prochlorococcus*, *Synechococcus*, a picoeukaryote, a small diatom,
- the code to generate histograms of measured count concentrations of *Prochlorococcus*, *Synechococcus*, their viruses and percentage of infected cells (Data from Carlson *et al.* 2022, https://www.nature.com/articles/s41564-022-01088-x).


All simulations were run on the Zaratan HPC cluster of the University of Maryland (https://hpcc.umd.edu/hpcc/zaratan.html)

Python libraries: numpy, matplotlib, math, sys, copy, scipy, random  
Python version: 3.12.9

R packages: scales  
R version 4.3.2

# 1. FILES DESCRIPTION

## 1.1. Python scripts:
For each file, pdf outputs are specified and *\*suffix\** refers to a suffix to identify parameters input of the simulation (see inside .py file and .sbatch files for inputs needed for each files). When specified, pdf outputs can depend on the input parameters: 
  - `generic_functions.py`: files containing functions shared by different scripts 
  - `SVZ_molar_model_phytotypes.py` : code to run SVZ simulations and generate coexistence diagrams
      - outputs:
        - SVZ_model_phi_latent_period_time_series_*\*suffix\**.pdf: time series for different adsorption rates
        - SVZ_model_phi_latent_period_*\*suffix\**.pdf: coexistence diagram and other features in the parameter space
  - `SVZ_functions.py` : functions specific to the SVZ model
  - `SVRZ_molar_model_phytotypes.py` : code to run SVRZ simulations and generate coexistence diagrams
    - outputs:
      - if dz2==0  (dz2 is an input parameter), SVRZ_model_phi_latent_period_time_series_coex_*\*suffix\**.pdf: if the quadratic mortality term of the zooplankton is 0: time series in the case of coexistence
      - SVRZ_model_phi_latent_period_time_series_*\*suffix\**.pdf: time series for different adsorption rates with phi/phir=10
      - SVRZ_model_phi_latent_period_time_series_full-res_*\*suffix\**.pdf: time series for different adsorption rates with phir=10 (full resistance)
      - SVRZ_model_phi_versus_phir_*\*suffix\**.pdf: coexistence diagram and other features in the parameter space
  - `SVRZ_functions.py` : functions specific to the SVRZ model
  - `SIVZ_molar_model_phytotypes.py` : code to run SIVZ simulations and generate coexistence diagrams
    - outputs:
      - if dz2==dv2==0 (dz2 and dv2 are input parameters), SIVZ_model_phi_latent_period_time_series_coex_*\*suffix\**.pdf: if the quadratic mortality terms are 0, time series in the case of coexistence
      - SIVZ_model_phi_latent_period_time_series_*\*suffix\**.pdf: time series for different adsorption rates, fixed latent period (life history trait model)
      - SIVZ_model_phi_latent_period_*\*suffix\**.pdf: coexistence diagram and other features in the parameter space
  - `SIVZ_functions.py` : functions specific to the SIVZ model
  - `SIVRZ_molar_model_phytotypes.py` : code to run SIVRZ simulations and generate coexistence diagrams
    - outputs:
      - SIVRZ_model_phi_latent_period_time_series_coex_*\*suffix\**.pdf: if the quadratic mortality terms are 0, time series in the case of coexistence
      - if param=='phir' or param=='epsr' (param is an input parameter), SIVRZ_model_phi_latent_period_time_series_*\*suffix\**.pdf, in case the resistance strength parameter space is explored, time series for different adsorption rates with phi/phir=10
      - if param=='phir' or param=='epsr', SIVRZ_model_phi_latent_period_time_series_full-res_*\*suffix\**.pdf', in case the resistance strength parameter space is explored, time series for different adsorption rates with phir=0 (full resistance)
      - Main pdf:
        - if param=='phir', SIVRZ_model_phi_versus_phir_*\*suffix\**.pdf, in the case extracellular resistance space is explored, coexistence diagram and other features in the parameter space
        - if param=='epsr', SIVRZ_model_eps_versus_epsr_*\*suffix\**.pdf, in the case extracellular resistance space is explored, coexistence diagram and other features in the parameter space
        - if param=='lp_phir' or param=='lp_epsr', in the case latent period space is explored (fixed intra or extra cellular resistance), coexistence diagram and other features in the parameter space
  Note: the SVZ_molar_model_phytotypes.py, SVRZ_molar_model_phytotypes.py SIVZ_molar_model_phytotypes.py, and SIVRZ_molar_model_phytotypes.py also save as .txt files some of the matrices they generate
  - `SIVRZ_functions.py` : functions specific to the SIVRZ model
  - `SIVZ_MCT.py`: code to run the SIVZ Modern coexistence theory analysis
    - outputs:
      - SIVZ_model_phi_latent_period_*\*suffix\**_MCT.pdf, effects from Modern Coexistence Theory (MCT, Ellner *et al.* 2019, https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13159) across the parameter space
  - `make_figures_scaled.py`: code to create figures with same scales (from different models): => figure 5, s6, s14, and s16
    - outputs:
      - Scaled_PVZ_main_models.pdf
      - Scaled_PVZ_FFT_and_Perdiod.pdf
      - Scaled_SIVZ_coexistence.pdf
      - Scaled_SIVZ_distances.pdf
  - `MCT_replicates_analysis.py`: generate figure3a-f (MCT analysis)
  - `SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py`: file to run the optimization of parameterization of the SIVZ and SIVRZ models with respect to target concentrations (grid search)
  - `analysis_optimization.py`: analyze the results of the optimization =>  generates data for table 1 and 2

## 1.2. R script
  - `histogram_abundances.R`: generate the histogram of count concentrations of *Prochlorococcus*, *Synechococcus*, their viruses and percentage of infected cells (Data from Carlson *et al.* 2022).
    - output: hist_abundances_Syn_Pro_virus_percentage_infected.pdf    

## 1.3. .sbatch files:
Each slurm sbatch file will need to be edited to conform to your HPC system. Note: .sbatch files will need module load for R on some HPC.
Each file contains a description of the parameters taken as inputs: .sbatch files are used to run simulations on a HPC cluster with slurm. 
  - `run_SVZ_molar_model_phytotypes.sbatch`: runs `SVZ_molar_model_phytotypes.py`
  - `run_SVRZ_molar_model_phytotypes.sbatch`: runs `SVRZ_molar_model_phytotypes.py`
  - `run_SIVZ_molar_model_phytotypes.sbatch`: runs `SIVZ_molar_model_phytotypes.py`
  - `run_SIVRZ_molar_model_phytotypes.sbatch`: runs `SIVRZ_molar_model_phytotypes.py`
  - `run_SIVZ_MCT.sbatch` runs `SIVZ_MCT.py`
  - `run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch`: runs `SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.py`
  - `run_analysis_optimization.sbatch`: runs `analysis_optimization.py`

## 1.4. .sh scripts: 
  For all .sh files, first add executable permissions: chmod +x name_of_file.sh
  - `run_coexistence_simulations.sh` : runs all simulations necessary to generate figure 2, 4, 6b-j, S5, S7, s8, s9, s10, s11, s12, s13, s17 and data for: figure 5, 3d,e, S6, s14, s16
  - `run_simulations_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sh`: used to run in parallel  the optimization of parameters, it divides the grid search in 50 equal chunks of parameters combinations e.g. if 10^6 parameter combination are tested, one job will test 10^6/50= 20000 combinations 
  - `concatenate_results_optimisation_params.sh`: used to concatenate result files from the optimization
  - for all .sh files, give permission to execute using chmod +x file.sh

## 1.5. .txt files:
  - `carlson_data_cyanobacteria.txt`: count concentration of *Prochlorococcus* and *Synechococcus* from Carlson *et al.* 2022
  - `carlson_data_virus_and_percentage_infected.txt`: count concentrations of cyanophages and percentage of infected cells (*Prochlorococcus* and *Synechococcus*) from Carlson *et al.* 2022
  - `coexistence_analysis_V_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic_Fluctuation_free_growth_rate.txt` and `coexistence_analysis_V_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic_Relative_non_linearity_I.txt` (necessary for example data to run `python make_figures_scaled.py main_models`)

# 2. RUN THE SIMULATIONS

In a HPC cluster (recommended), or locally (probably ~6-7 days of computation):  

## 2.1. Run simulations of dynamics across the parameter space:  
WARNING: this step will overwrite example data stored in `model_data/` and the files `coexistence_analysis_V_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic_Fluctuation_free_growth_rate.txt` and `coexistence_analysis_V_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic_Relative_non_linearity_I.txt`  
&nbsp;&nbsp;&nbsp;&nbsp;`module load python`  
&nbsp;&nbsp;&nbsp;&nbsp;- run scripts using slurm: `./run_coexistence_simulations.sh` OR  
&nbsp;&nbsp;&nbsp;&nbsp;- run scripts using python in your terminal: `run_coexistence_simulations_python_commands.sh` (~5.5 days of computing)  
&nbsp;&nbsp;&nbsp;&nbsp; This step will run all simulations necessary to generate figures of the coexistence diagram and some data that needs to be rescaled to generate some figures (already preloaded in `model_data/`).

&nbsp;&nbsp;&nbsp;&nbsp;**Run times:**  
&nbsp;&nbsp;&nbsp;&nbsp;- SVZ and SVRZ simulations across the adsorption rate parameter space: ~3 minutes  
&nbsp;&nbsp;&nbsp;&nbsp;- SIVZ and SIVRZ simulations across adsorption rate and latent period: ~2h30 to 3 hours  
&nbsp;&nbsp;&nbsp;&nbsp;- SVRZ and SIVRZ across adsorption rate and resistance strength: ~3h30 to 5 hours  
&nbsp;&nbsp;&nbsp;&nbsp;- MCT jobs: 6 hours each

&nbsp;&nbsp;&nbsp;&nbsp;Before running all together, consider checking each job (SVZ, SIVZ, SVRZ, and SIVRZ) individually to ensure it works.

## 2.2. Generate MCT figures (figure 3a-f):
&nbsp;&nbsp;&nbsp;&nbsp;`python MCT_replicates_analysis.py`  
&nbsp;&nbsp;&nbsp;&nbsp; This will analyse the replicates of the MCT analysis and output figure 3a-f and the sensitivity og the Chesson criterion

## 2.3. Generate scaled figures:  
&nbsp;&nbsp;&nbsp;&nbsp;`python make_figures_scaled.py main_models`  
&nbsp;&nbsp;&nbsp;&nbsp;This will generate figures that need to be scaled (necessary data previously saved in the `model_data/` folder).

## 2.4. Generate figures of distributions of measured abundances of *Prochlorococcus* and *Synechococcus*, their virus and the percentage of infected cells from Carlson *et al.* 2022.
&nbsp;&nbsp;&nbsp;&nbsp;`Rscript histogram_abundances.R`

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
Note that the message `cat: results_optimization_params/SIVZ_intracellular_res_optimization_*organism*_*num*.txt: No such file or directory` is normal (no realistic concentrations found in this chunk)

## 2.8. Run the analysis of the optimization:  
&nbsp;&nbsp;&nbsp;&nbsp;`sbatch run_analysis_optimization.sbatch 0`  
&nbsp;&nbsp;&nbsp;&nbsp; This will extract the best and 200 best parameter combinations for each phytoplankton type

# 3. MANUSCRIPT FIGURE PANELS AND TABLE DATA (PDFs and TXTs) 

After running all simulations, the following .pdf and .txt files contain the figure panels and data to generate figure 2 to 6, figure S5 to S17 and Table 1 and 2.  

## 3.1. Figures
Pages of figure panels are reported in this section. Figures were assembled using Inkscape 1.4.2
- **Figure 2**
  - 2a: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.pdf
  - 2b: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - 2c: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_no-dz2_m2-1800.0_mesotrophic.pdf
  - 2d: page 1 of SVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 2e: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.pdf
  - 2f: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - 2g: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_no-dz2_mesotrophic.pdf
  - 2h: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 2i: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf
  - 2j: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir-5.0.pdf
  - 2k: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_phir-5.0.pdf
  - 2l: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir-5.0.pdf
  - 2m: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_lp_phi-ratio-5.0_mesotrophic.pdf
  - 2n: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - 2o: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 2p: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 2q: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-0.0.pdf
  - 2r: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir-0.0.pdf
  - 2s: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_phir-0.0.pdf
  - 2t: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir-0.0.pdf
  - 2u: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_lp_phi-ratio-0.0_mesotrophic.pdf
  - 2v: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-0.0_mesotrophic.pdf
  - 2w: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf
  - 2x: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf
- **Figure 3**
  - 3a: page 20 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_no-dz2_no-dv2_mesotrophic_MCT.pdf
  - 3b: page 19 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_no-dz2_no-dv2_mesotrophic_MCT.pdf
  - 3c: page 15 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_no-dz2_no-dv2_mesotrophic_MCT.pdf
  - 3d: page 3 of Scaled_SIVZ_coexistence.pdf
  - 3e: page 4 of Scaled_SIVZ_coexistence.pdf
  - black MASK: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.pdf
- **Figure 4**
  - 4a: page 6 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - 4b: page 10 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - 4c: page 14 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - 4d: page 18 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - 4e: page 6 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - 4f: page 10 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - 4g: page 14 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - 4h: page 18 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - 4i: page 6 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - 4j: page 10 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - 4k: page 14 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - 4l: page 18 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - 4m: page 6 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - 4n: page 10 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - 4o: page 14 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - 4p: page 18 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
- **Figure 5**
  - 5a: page 1 of Scaled_PVZ_main_models.pdf
  - 5b: page 2 of Scaled_PVZ_main_models.pdf
  - 5c: page 3 of Scaled_PVZ_main_models.pdf
  - 5d: page 4 of Scaled_PVZ_main_models.pdf
  - 5e: page 5 of Scaled_PVZ_main_models.pdf
  - 5f: page 6 of Scaled_PVZ_main_models.pdf
  - 5g: page 7 of Scaled_PVZ_main_models.pdf
  - 5h: page 8 of Scaled_PVZ_main_models.pdf
  - 5i: page 9 of Scaled_PVZ_main_models.pdf
  - 5j: page 19 of Scaled_PVZ_main_models.pdf
  - 5k: page 20 of Scaled_PVZ_main_models.pdf
  - 5l: page 21 of Scaled_PVZ_main_models.pdf
  - 5m: page 10 of Scaled_PVZ_main_models.pdf
  - 5n: page 11 of Scaled_PVZ_main_models.pdf
  - 5o: page 12 of Scaled_PVZ_main_models.pdf
  - 5p: page 34 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 5q: page 43 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 5r: page 43 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf
  - 5s: page 49 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 5t: page 10 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 5u: page 10 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
- **Figure 6**
  - 6b: page 47 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 6c: page 60 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 6d: page 60 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 6e: page 51 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 6f: page 64 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 6g: page 64 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 6h: page 56 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - 6i: page 69 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - 6j: page 69 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf

## 3.2. Tables

- **Table 1**
  - Best parameters: SIVZ_optimum_fit_*\*phyto-type\**_intracellular.txt where *\*phyto-type\** is Prochlorococcus, Synechococcus, Eukaryote and Diatom
  - Parameter ranges: SIVZ_range_fit_*\*phyto-type\**_intracellular.txt where *\*phyto-type\** is Prochlorococcus, Synechococcus, Eukaryote and Diatom
- **Table 2**
  - Minimum total error: SIVZ_optimum_fit_*\*phyto-type\**_intracellular.txt
  - Error range: SIVZ_range_fit_*\*phyto-type\**_intracellular.txt
  - Reached concentrations: SIVZ_optimum_fit_errors_tracers_*\*phyto-type\**_intracellular.txt

## 3.3. Supplementary figures

- **Figure S1** to **Figure S4**: see trait_data folder
- **Figure S5** SIVZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.pdf 
  - S5a: pages 17, 19 and 20
  - S5b: pages 48, 50 and 51
  - WARNING: this figure might differ from one system to another (numerical instability)
- **Figure S6** Scaled_PVZ_FFT_and_Period.pdf
  - S6a to S6h: pages 1 to 8
- **Figure S7**
  - S7a: pages 1 and 3 of SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf
  - S7b: pages 17 and 19 of SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf
  - S7c: pages 1 and 3 of SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-fullres.pdf
  - S7d: pages 17 and 19 of SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-fullres.pdf
- **Figure S8**
  - S8a: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir.pdf
  - S8b: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir.pdf
  - S8c: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800_no-dz2_mesotrophic_phir.pdf
  - S8d: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800_mesotrophic_phir.pdf
  - S8e: page 1 of SIVRZ_model_phi_versus_phir_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_phir_mesotrophic.pdf
  - S8f: page 1 of SIVRZ_model_phi_versus_phir_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_phir_mesotrophic.pdf
  - S8g: page 1 of SIVRZ_model_phi_versus_phir_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_no-dz2_phir_mesotrophic.pdf
  - S8h: page 1 of SIVRZ_model_phi_versus_phir_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S8i: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_epsr.pdf
  - S8j: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_epsr.pdf
  - S8k: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800_no-dz2_mesotrophic_epsr.pdf
  - S8l: page 1 of SVRZ_model_phi_versus_phir_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800_mesotrophic_epsr.pdf
  - S8m: page 1 of SIVRZ_model_eps_versus_epsr_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_phir_mesotrophic.pdf
  - S8n: page 1 of SIVRZ_model_eps_versus_epsr_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_phir_mesotrophic.pdf
  - S8o: page 1 of SIVRZ_model_eps_versus_epsr_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_no-dz2_phir_mesotrophic.pdf
  - S8p: page 1 of SIVRZ_model_eps_versus_epsr_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
- **Figure S9**
  - S9a: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - S9b: page 7 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - S9c: page 1 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S9d: page 7 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S9e: page 1 of SIVZ_model_phi_latent_period_Synechococcus_BS15.0_LOI0_mesotrophic.pdf
  - S9f: page 7 of SIVZ_model_phi_latent_period_Synechococcus_BS15.0_LOI0_mesotrophic.pdf
  - S9g: page 1 of SIVZ_model_phi_latent_period_Synechococcus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S9h: page 7 of SIVZ_model_phi_latent_period_Synechococcus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S9i: page 1 of SIVZ_model_phi_latent_period_Eukaryote_BS15.0_LOI0_mesotrophic.pdf
  - S9j: page 7 of SIVZ_model_phi_latent_period_Eukaryote_BS15.0_LOI0_mesotrophic.pdf
  - S9k: page 1 of SIVZ_model_phi_latent_period_Eukaryote_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S9l: page 7 of SIVZ_model_phi_latent_period_Eukaryote_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S9m: page 1 of SIVZ_model_phi_latent_period_Diatom_BS15.0_LOI0_mesotrophic.pdf
  - S9n: page 7 of SIVZ_model_phi_latent_period_Diatom_BS15.0_LOI0_mesotrophic.pdf
  - S9o: page 1 of SIVZ_model_phi_latent_period_Diatom_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S9p: page 7 of SIVZ_model_phi_latent_period_Diatom_BS15.0_LOI0_m2-1800_mesotrophic.pdf
- **Figure S10**
  - S10a to S10e: pages 2 to 6 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - S10f to S10j: pages 2 ro 6 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
- **Figure S11**
  - S9a: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9b: page 13 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9c: page 1 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9d: page 13 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9e: page 1 of SIVRZ_model_phi_versus_latent_period_Synechococcus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9f: page 13 of SIVRZ_model_phi_versus_latent_period_Synechococcus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9g: page 1 of SIVRZ_model_phi_versus_latent_period_Synechococcus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9h: page 13 of SIVRZ_model_phi_versus_latent_period_Synechococcus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9i: page 1 of SIVRZ_model_phi_versus_latent_period_Eukaryote_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9j: page 13 of SIVRZ_model_phi_versus_latent_period_Eukaryote_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9k: page 1 of SIVRZ_model_phi_versus_latent_period_Eukaryote_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9l: page 13 of SIVRZ_model_phi_versus_latent_period_Eukaryote_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9m: page 1 of SIVRZ_model_phi_versus_latent_period_Diatom_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9n: page 13 of SIVRZ_model_phi_versus_latent_period_Diatom_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9o: page 1 of SIVRZ_model_phi_versus_latent_period_Diatom_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S9p: page 13 of SIVRZ_model_phi_versus_latent_period_Diatom_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
- **Figure S12**
  - S12a to S12k: pages 2 to 12 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S12l to S12v: pages 2 to 12 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
- **Figure S13**
  - S13a: page 5 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - S13b: page 9 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - S13c: page 13 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - S13d: page 17 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf
  - S13e: page 5 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S13f: page 9 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S13g: page 13 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S13h: page 17 of SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf
  - S13i: page 5 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S13j: page 9 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S13k: page 13 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S13l: page 17 of SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S13m: page 5 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S13n: page 9 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S13o: page 13 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
  - S13p: page 17 of SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf
- **Figure S14**
  - S14a to S14f: pages 13 to 18 of Scaled_PVZ_main_models.pdf
  - S14g and S14h: pages 22 and 23 of Scaled_PVZ_main_models.pdf
  - S14i: page 44 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S14j: page 44 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf
- **Figure S15**
  - S15a to S15d: pages 1,2,4,5 of hist_abundances_Syn_Pro_virus_percentage_infected.pdf
  - S15e: page 3 of hist_abundances_Syn_Pro_virus_percentage_infected.pdf
- **Figure S16**
  - S16a to S16i: pages 1 to 9 of Scaled_SIVZ_distances.pdf
- **Figure S17**
  - S17a: page 48 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - S17b: page 61 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S17c: page 61 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S17d: page 49 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - S17e: page 62 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S17f: page 62 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S17g: page 50 of SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf
  - S17h: page 63 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
  - S17i: page 63 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf

