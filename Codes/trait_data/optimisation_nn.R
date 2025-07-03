best_models_nn<-function(id,  df_full){
  # model training function
  optimisation_nn <- function(i, nnetGrid, variables, df_full, id){
    flag <- TRUE
    nn_model <- NULL
    # Fitting the model on the whole dataset
    set.seed(nnetGrid[i,4])
    possibleError <- tryCatch(
      nn_model <- nnet::nnet(y=df_full[,id], x=df_full[,variables],
                             data = df_full, size=nnetGrid[i,1], decay=nnetGrid[i,2],
                             maxit=nnetGrid[i,3], trace=F, linout=T),
      error=function(e) flag <- FALSE
    )

    # cross validation function: fit the model on a subset and compute evaluation metrics on the test set
    cv_nn <- function(sample,i, id){
      df3 <- df_full[sample,]
      test_set <- !(c(1:nrow(df_full)) %in% sample)
      set.seed(nnetGrid[i,4])
      nn_model_cv <- nnet::nnet(y=df3[,id], x=df3[,variables],
                                data = df3, size=nnetGrid[i,1],
                                decay=nnetGrid[i,2], maxit=nnetGrid[i,3],
                                trace=F, linout=T)
      preds <- stats::predict(nn_model_cv, df_full[test_set,variables], type='raw')
      rms <- Metrics::rmse(actual=df_full[test_set,id], predicted=preds)
      co <- cor(df_full[test_set,id],preds, method = 'pearson' )
      return(c(rms, co))
    }

    # divide the data set LOOCV
    set.seed(42)
    to_test=1:dim(df_full)[1]
    to_test=to_test[!(to_test %in% to_rem)]
    samples_list <- rep(list(NULL), length(to_test))
    co=1
    for (u in to_test){
      samples <- 1:dim(df_full)[1]
      samples=samples[samples!=u]
      samples_list[[co]] <- samples
      co=co+1
    }
    # run LOOCV
    score_list <- lapply(samples_list, FUN = cv_nn, i=i, id=id)

     # collect model performance
    rmse_list <- NULL
    cor_list <- NULL
    for (vector in score_list){
      rmse_list <- append(rmse_list,vector[1])
      cor_list <- append(cor_list,vector[2])
    }
    RMSEs <- mean(rmse_list, na.rm=T)
    CORs <- mean(cor_list, na.rm=T)
    if (flag == FALSE || is.null(nn_model)) {
      RMSEs <- NA
      CORs <- NA
    }
    optimisation_nn <- list(RMSEs, CORs, nn_model)
  }

  # main loop over the hyperparameter space
  m_core <- NULL
  m_rmse <- NULL
  mods <- NULL
  for (i in 1:dim(nnetGrid)[1]){ #
    set.seed(42)
    # LOOCV for a given combination of hyperparameter
    vec <- optimisation_nn(i, nnetGrid = nnetGrid , variables=variables, df_full = df_full, id=id)
    m_rmse <- append(m_rmse, vec[[1]])
    m_core <- append(m_core, vec[[2]])
    mods[[i]] <- vec[[3]]
  }
  # Finding the parameters combination for which the mean auc of the cross-validation is maximized
  g <- which.min(m_rmse)
  if (length(g)>0){
    best_model <- c(nnetGrid[g,1], nnetGrid[g,2],nnetGrid[g,3], nnetGrid[g,4], m_rmse[g], m_core[g])
  } else{
    best_model <- c(NA,NA,NA,NA, NA, NA)
  }

  b_mod <- mods[[g]]
  # return the best model and its performance
  return(list(best_model, b_mod))
}
  
