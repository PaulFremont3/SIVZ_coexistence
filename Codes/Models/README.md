In this directory: all python and R codes to:
- run simulations of the SVZ, SIVZ, SVRZ and SIVRZ models,
- example data in the folder `model_data/` from Figure 4, Figure 7a-o, Figure S7, Figure s9, Figure s10, Figure s12, Figure s14, and Figure s16 are stored in the folder `model_data` so these figures can be directly reproduced by running:  
   - `python make_figures_scaled.py main_models` for Figure 7a-o, Figure S7, Figure s14 and Figure s16: generating the files `Scaled_PVZ_FFT_and_Period.pdf`,  `Scaled_PVZ_main_models.pdf`,  `Scaled_SIVZ_coexistence.pdf`, and  `Scaled_SIVZ_distances.pdf` in <10 seconds)
   -  `python make_state_diagrams_figures.py` for Figure 4, Figure s9, Figure s12 and Figure s14: generating the files `State_diagrams_main_models.pdf`, `State_diagrams_resistance_models.pdf`, `State_diagrams_phytoplankton_types_models.pdf`, `State_diagrams_phytoplankton_types_resistance_models.pdf`, `State_diagrams_phytoplankton_types_th.pdf`, and `State_diagrams_phytoplankton_types_resistance_th.pdf`)  
- to reproduce other figures, simulations have to be run
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
  - `make_figures_scaled.py`: code to create figures with same scales (from different models): => figure 7a-o, s7, s14, and s16
    - outputs:
      - `Scaled_PVZ_main_models.pdf`
      - `Scaled_PVZ_FFT_and_Perdiod.pdf`
      - `Scaled_SIVZ_coexistence.pdf`
      - `Scaled_SIVZ_distances.pdf`
  - `make_state_diagrams_figures.py`: code to create state diagram figures: Figure 4, Figure s9, Figure s12 and Figure s14
     - outputs:
       - `State_diagrams_main_models.pdf`
       - `State_diagrams_resistance_models.pdf`
       - `State_diagrams_phytoplankton_types_models.pdf`
       - `State_diagrams_phytoplankton_types_resistance_models.pdf`
       - `State_diagrams_phytoplankton_types_th.pdf`
       - `State_diagrams_phytoplankton_types_resistance_th.pdf`
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
  - `run_coexistence_simulations.sh` : runs all simulations necessary to generate figures 3, 5, 6, 7p-u, 8, s5, s6, s8, s11, s13, s15, s16i-j, s18 and data for 4, 7a-o, s7, s9, s10, s12, s14, s16a-h
  - `run_coexistence_simulations_python_commands.sh`: runs all simulations necessary to generate figures
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
&nbsp;&nbsp;&nbsp;&nbsp; This will analyse the replicates of the MCT analysis and output pdfs for figure 5 and the sensitivity of the Chesson criterion

## 2.3. Generate scaled figures:  
&nbsp;&nbsp;&nbsp;&nbsp;`python make_figures_scaled.py main_models`  
&nbsp;&nbsp;&nbsp;&nbsp;This will generate figures that need to be scaled (necessary data previously saved in the `model_data/` folder, running will override).

## 2.3. Generate state diagram figures:  
&nbsp;&nbsp;&nbsp;&nbsp;`python make_state_diagrams_figures.py main_models`  
&nbsp;&nbsp;&nbsp;&nbsp;This will generate all state diagram figures (necessary data previously saved in the `model_data/` folder, running will override).

## 2.5. Generate figures of distributions of measured abundances of *Prochlorococcus* and *Synechococcus*, their virus and the percentage of infected cells from Carlson *et al.* 2022.
&nbsp;&nbsp;&nbsp;&nbsp;`Rscript histogram_abundances.R`

## 2.6. Create a subfolder to store parameter optimization results:  
&nbsp;&nbsp;&nbsp;&nbsp;`mkdir results_optimization_params/`

## 2.7. Grid search to optimize parameters to target concentrations:  
&nbsp;&nbsp;&nbsp;&nbsp;`./run_simulations_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sh 0`  
&nbsp;&nbsp;&nbsp;&nbsp;This will run optimization of parameters to minimize distance to target concentrations for the four phytoplankton types considered in the study:  
&nbsp;&nbsp;&nbsp;&nbsp;a small diatom, a picoeukaryote, a *Synechococcus*, and a *Prochlorococcus*.  
&nbsp;&nbsp;&nbsp;&nbsp;Results are stored in the `results_optimization_params/` folder.

## 2.8. Concatenate files of optimization results:  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Synechococcus 0`  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Prochlorochoccus 0`  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Eukaryote 0`  
&nbsp;&nbsp;&nbsp;&nbsp;`./concatenate_results_optimisation_params.sh SIVZ intracellular Diatom 0`
Note that the message `cat: results_optimization_params/SIVZ_intracellular_res_optimization_*organism*_*num*.txt: No such file or directory` is normal (no realistic concentrations found in this chunk)

## 2.9. Run the analysis of the optimization:  
&nbsp;&nbsp;&nbsp;&nbsp;`sbatch run_analysis_optimization.sbatch 0`  
&nbsp;&nbsp;&nbsp;&nbsp; This will extract the best and 200 best parameter combinations for each phytoplankton type

# 3. MANUSCRIPT FIGURE PANELS AND TABLE DATA (PDFs and TXTs) 

After running all simulations, the following .pdf and .txt files contain the figure panels and data to generate figure 2 to 6, figure S5 to S17 and Table 1 and 2.  

## 3.1. Figures
Pages of figure panels are reported in this section. Figures were assembled using Inkscape 1.4.2
- **Figure 3**
   -  3a-d: page 6, 10, 14 and 18 of `SVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.pdf`
   -  3e,g,h: page 6, 10, 18 of `SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic`
   -  3f: page 24 of `SIVZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic` (WARNIGN: this may vary from one system to another due to numerical instabilities)
   -  3i,k,l: age 6, 10, 18 of `SVRZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf`
   -  3j: page 10 of `SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf`
- **Figure 4**: a to x: page 1 to 24 of `State_diagrams_main_models.pdf`
- **Figure 5**: a to f: page 2,3,5,11,7,9 of `SIVZ_model_phi_latent_period_replicates_consensus_MCT.pdf`
- **Figure 6**
  - 6a-d: page 6, 10, 14, 18 of `SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf`
  - 6e-h: page 6, 10, 14, 18 of `SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf`
  - 6i-l: page 6, 10, 14, 18 of `SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf`
  - 6m-p: page 6, 10, 14, 18 of `SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf`
- **Figure 7**
  - 7a-o: page 1-9, 19-21, 11, 12 of `Scaled_PVZ_main_models.pdf`
  - 7p: page 34 of `SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf`
  - 7q: page 43 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`
  - 7r: page 43 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf`
  - 7s: page 49 of `SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf`
  - 7t: page 10 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`
  - 7u: page 10 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`
- **Figure 8**
  - 6b-d: page 1 to 3 of `Scaled_SIVZ_distances_tot.pdf`
  - 6e,h: page 51 and 56 of `SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf`
  - 6f,i: page 64 and 69 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`
  - 6g,j: page 64 and 69 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`

## 3.2. Tables

- **Table 1**
  - Best parameters: SIVZ_optimum_fit_*\*phyto-type\**_intracellular.txt where *\*phyto-type\** is Prochlorococcus, Synechococcus, Eukaryote and Diatom
  - Parameter ranges: SIVZ_range_fit_*\*phyto-type\**_intracellular.txt where *\*phyto-type\** is Prochlorococcus, Synechococcus, Eukaryote and Diatom
- **Table 2**
  - Minimum total error: SIVZ_optimum_fit_*\*phyto-type\**_intracellular.txt
  - Error range: SIVZ_range_fit_*\*phyto-type\**_intracellular.txt
  - Reached concentrations: SIVZ_optimum_fit_errors_tracers_*\*phyto-type\**_intracellular.txt

## 3.3. Supplementary figures

- **Figure S1** to **Figure S4**: see `trait_data/` folder
- **Figure S5**
   -  s5a-d: page 5, 9, 13 and 17 of `SVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.pdf`
   -  s5e,g,h: page 5, 9, 17 of `SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic`
   -  s5f: page 25 of `SIVZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic` (WARNING: this may vary from one system to another due to numerical instabilities)
   -  s5i,k,l: page 5, 9, 17 of `SVRZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf`
   -  s5j: page 11 of `SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf`
- **Figure S6**: a and b: page 45, 47, 50, 244, 246 and 249 of `SIVZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.pdf` (WARNING: this may vary from one system to another due to numerical instabilities)
- **Figure S7**: page 3 and 4 of `Scaled_SIVZ_coexistence.pdf`
- **Figure S8**
   - s8a and b: page 1,3, 31 and 33 of `SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.pdf`
   - s8c and d: page 1,3, 31 and 33 of `SVRZ_model_phi_latent_period_time_series_coex_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-fullres.pdf`
- **Figure S9**: s9a to p: page 1 to 16 of `State_diagrams_resistance_models.pdf`
- **Figure S10**:
   - s10a,c,e,g,i,k,m,o: page 1 to 8 of `State_diagrams_phytoplankton_types_models.pdf`
   - s10b,d,f,h,j,l,n,p: page 1 to 8 of `State_diagrams_phytoplankton_types_models_th.pdf`
- **Figure S11**:
   - S11a to S11e: pages 2 to 6 of `SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf`
   - S11f to S11j: pages 2 ro 6 of `SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf`
- **Figure S12**:
   - s12a,c,e,g,i,k,m,o: page 1 to 8 of `State_diagrams_phytoplankton_types_resistance_models.pdf`
   - s12b,d,f,h,j,l,n,p: page 1 to 8 of `State_diagrams_phytoplankton_types_resistance_th.pdf`
- **Figure S13**
  - S13a to S13k: pages 2 to 12 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.pdf
  - S13l to S13v: pages 2 to 12 of SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf
- **Figure S14** Scaled_PVZ_FFT_and_Period.pdf
  - S6a to S6h: pages 1 to 8
- **Figure S15**
  - s15a-d: page 5, 9, 13, 17 of `SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_mesotrophic.pdf`
  - s15e-h: page 5, 9, 13, 17 of `SIVZ_model_phi_latent_period_time_series_Prochlorochoccus_BS15.0_LOI0_m2-1800_mesotrophic.pdf`
  - s15i-l: page 5, 9, 13, 17 of `SIVRZ_model_phi_latent_period_time_series_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf`
  - s15m-p: page 5, 9, 13, 17 of `SIVRZ_model_phi_latent_period_time_series_full-res_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_phir_mesotrophic.pdf`
- **Figure S16**
  - S16a to S16f: pages 13 to 18 of `Scaled_PVZ_main_models.pdf`
  - S16g and S16h: pages 22 and 23 of `Scaled_PVZ_main_models.pdf`
  - S16i: page 44 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`
  - S16j: page 44 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-0.0_mesotrophic.pdf`
- **Figure S17**: see `trait_data/` folder
- **Figure S18**: s18a-i: page 1 to 9 of `Scaled_SIVZ_distances.pdf`
- **Figure S19**
  - S17a,d,g: page 48,49,63 of `SIVZ_model_phi_latent_period_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.pdf`
  - S17b,e,h: page 61,62,63 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`
  - S17c,f,i: age 61,62,63 of `SIVRZ_model_phi_versus_latent_period_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800_lp_phi-ratio-5.0_mesotrophic.pdf`
 
