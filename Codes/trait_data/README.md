This folder contains all R codes necessary to generate life history traits models.

- R scripts
    - optimisation_lm.R: cross validation function for the linear model
    - optimisation_nn.R: cross validation function for the single layer neural network model
    - optimisation_rf.R: cross validation function for the random forest model
    - optimisation_gam.R: cross validation function for the gam model
    - best_model_ml.R: generic code to launch either the optimisation of lm, nn, rf or gam model and generates figure s1 to s4
    - traits_for_simulation.R save final predictions of the model
 
- .sbatch file
  - run_best_model_ml.sbatch: to run the best_model_ml.R script (see arguments inside)
