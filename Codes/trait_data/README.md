This folder contains all R codes necessary to generate life history traits models. Note that all data is already in the file so you can directly run models in the `Models/` folder.

All simulations were run on the Zaratan HPC cluster of the University of Maryland (https://hpcc.umd.edu/hpcc/zaratan.html).

# 1. FILES DESCRIPTION

## 1.1. R scripts
  - optimisation_lm.R: cross validation function for the linear model
  - optimisation_nn.R: cross validation function for the single layer neural network model
  - optimisation_rf.R: cross validation function for the random forest model
  - optimisation_gam.R: cross validation function for the gam model
  - best_model_ml.R: generic code to launch either the optimisation of lm, nn, rf or gam model and generates figure s1 to s4
  - traits_for_simulation.R save final predictions of the model
 
## 1.2 .sbatch file
  - run_best_model_ml.sbatch: to run the best_model_ml.R script (see arguments inside)

## 1.3. .txt files
  - Vs_5.txt: 400 random volumes of 100 diatom, 100 eukaryotes, 100 Synechococcus and 100 Prochlorococcus
  - mumax_dutkiewicz_5.txt: maximum growth rates of the 400 random phytoplankton based on allometric relationship Dutkiewicz *et al.* 2020
  - Nc_dutkiewicz_5.txt: half saturation constant of the 400 random phytoplankton based on allometric relationship from Dutkiewicz *et al.* 2020 (https://bg.copernicus.org/articles/17/609/2020/)
  - model_burst_size_nn-gam.txt: burst size model generated from the life history trait model of this study (nn for burst size of eukaryotes, Synechococcus and Prochlorococcus and gam for burst size of diatom)
  - model_latent_period_nn-gam.txt: burst size model generated from the life history trait model of this study (nn for latent period)
Practically in the study, only the life history trait of the smallest phytoplankton of each group is used.

# 2. RUN THE SIMULATIONS
## 2.1. Train modles of viral life history traits
- Burst size:  
`sbatch run_best_model_ml.sbatch nn BS`  
`sbatch run_best_model_ml.sbatch gam BS`  
`sbatch run_best_model_ml.sbatch lm BS`  
`sbatch run_best_model_ml.sbatch rf BS`  
- Latent period:
`sbatch run_best_model_ml.sbatch nn BS`  
`sbatch run_best_model_ml.sbatch gam BS`  
`sbatch run_best_model_ml.sbatch lm BS`  
`sbatch run_best_model_ml.sbatch rf BS`  
## 2.2. Generate life history traits of phytoplankton
Compile the c++ file generate_traits.cpp:  
`g++ generate_traits.cpp -o generate_traits.exe`  
Generate the trait data:
`./generate_traits.exe`

# 3. PDF FILES WITH THE FIGURES IN THE PAPER
