This folder contains all R codes necessary to generate life history traits models.

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

# 2. RUN THE SIMULATIONS

# 3. PDF FILES WITH THE FIGURES IN THE PAPER
