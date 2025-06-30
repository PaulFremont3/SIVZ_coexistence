This folder contains all R codes necessary to generate life history traits models. Note that all data is already in the file so you can directly run models in the `Models/` folder. When running codes in the `trait_data` folder, it will replace the original files but they should be identical.

All simulations were run on the Zaratan HPC cluster of the University of Maryland (https://hpcc.umd.edu/hpcc/zaratan.html).

R libraries: mgcv, nnet, Metrics
R version: 3.3.2

# 1. FILES DESCRIPTION

## 1.1. R scripts
  - `optimisation_lm.R`: cross validation function looped through the hyperparameter space for the linear model
  - `optimisation_nn.R`: cross validation function looped through the hyperparameter space for the single layer neural network model
  - `optimisation_rf.R`: cross validation function looped through the hyperparameter space for the random forest model
  - `optimisation_gam.R`: cross validation function looped through the hyperparameter space for the gam model
  - `best_model_ml.R`: generic code to launch either the optimisation of lm, nn, rf or gam model and generates figure s1 to s4
    - output:
      - pred_vs_data_*\*variable\**__*\*model\**.pdf:  model versus data pdf
      - model_vs_data_*\*variable\**_host_volume_*\*model\**.pdf: model and data versus volume of host pdf
  - `traits_for_simulation.R` save final predictions of the model
    - output:
      - trait_values_size_structured_foodweb.pdf: pdf with the trait values for the 400 phytoplankton types (see .cpp). In this study, only the smallest of each phytoplankton groups are used.

## 1.2. .cpp script
 - `generate_traits.cpp` (goes with `functions.h`): cpp file to generate randomly 400 phytoplankton types and their life history traits (based on their volumes): 100 diatom, 100 eukaryotes, 100 Synechococcus and 100 Prochlorococcus
 
## 1.3. .sbatch file
  - `run_best_model_ml.sbatch`: to run the best_model_ml.R script (see arguments inside)

## 1.4. .txt files
  - `Vs_5.txt`: 400 random volumes of 100 diatom, 100 eukaryotes, 100 Synechococcus and 100 Prochlorococcus
  - `mumax_dutkiewicz_5.txt`: maximum growth rates of the 400 random phytoplankton based on allometric relationship Dutkiewicz *et al.* 2020
  - `Nc_dutkiewicz_5.txt`: half saturation constant of the 400 random phytoplankton based on allometric relationship from Dutkiewicz *et al.* 2020 (https://bg.copernicus.org/articles/17/609/2020/)
  - `model_burst_size_nn-gam.txt`: burst size model generated from the life history trait model of this study (nn for burst size of eukaryotes, Synechococcus and Prochlorococcus and gam for burst size of diatom)
  - `model_latent_period_nn-gam.txt`: burst size model generated from the life history trait model of this study (nn for latent period)
Practically in the study, only the life history trait of the smallest phytoplankton of each group is used.
  - `edwards_2018.txt`: contains data from Edwards *et al.* 2018, necessary to train the models

# 2. RUN THE SIMULATIONS
## 2.1. Train models of viral life history traits
- Burst size:  
`sbatch run_best_model_ml.sbatch nn BS`  
`sbatch run_best_model_ml.sbatch gam BS`  
`sbatch run_best_model_ml.sbatch lm BS`  
`sbatch run_best_model_ml.sbatch rf BS`  
- Latent period:  
`sbatch run_best_model_ml.sbatch nn LP`  
`sbatch run_best_model_ml.sbatch gam LP`  
`sbatch run_best_model_ml.sbatch lm LP`  
`sbatch run_best_model_ml.sbatch rf LP`  
## 2.2. Generate life history traits data of phytoplankton
Compile the c++ script `generate_traits.cpp`:  
`g++ generate_traits.cpp -o generate_traits.exe`  
Generate the trait data:
`./generate_traits.exe 5` (5 is the random seed)
## 2.3. Generate life history traits data of viruses
`Rscript traits_for_simulation.R`

# 3. MANUSCRIPT FIGURE PANELS (PDFs)
- **Figure S1**
  - S1a: page 4 of model_vs_data_BS_host_volume_lm.pdf
  - S1b: page 4 of model_vs_data_BS_host_volume_gam.pdf
  - S1c: page 4 of model_vs_data_BS_host_volume_rf.pdf
  - S1d: page 4 of model_vs_data_BS_host_volume_nn.pdf
- **Figure S2**
  - S2a: page 1 of pred_vs_data_BS_lm.pdf
  - S2b: page 1 of pred_vs_data_BS_gam.pdf
  - S2c: page 1 of pred_vs_data_BS_rf.pdf
  - S2d: page 1 of pred_vs_data_BS_nn.pdf
- **Figure S3**
  - S3a: page 4 of model_vs_data_LP_host_volume_lm.pdf
  - S3b: page 4 of model_vs_data_LP_host_volume_gam.pdf
  - S3c: page 4 of model_vs_data_LP_host_volume_rf.pdf
  - S3d: page 4 of model_vs_data_LP_host_volume_nn.pdf
- **Figure S4**
  - S4a: page 1 of pred_vs_data_LP_lm.pdf
  - S4b: page 1 of pred_vs_data_LP_gam.pdf
  - S4c: page 1 of pred_vs_data_LP_rf.pdf
  - S4d: page 1 of pred_vs_data_LP_nn.pdf
