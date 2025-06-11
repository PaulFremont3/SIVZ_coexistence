best_models_rf<-function(id,  df_full){
  optimisation_rf <- function(i, rfGrid, variables, df_full, id){
    #print(ido)
    flag <- TRUE
    rf_model <- NULL
    # Fitting the model on the whole dataset
    set.seed(42)
    possibleError <- tryCatch(
      rf_model <- randomForest::randomForest(y=df_full[,id], x=df_full[,variables],data = df_full, importance=T,
                                             proximity =T, mtry=rfGrid[i,1], ntree = rfGrid[i,2]),
      error=function(e) flag <- FALSE
    )
    #test3<- nn_model$fitted.values[,2]
    #TSSs <- 0
    # Performing 15 cross-validations to calculate the mean AUC (area under ROC curve) which is the parameter we chose to optimize
    cv_rf <- function(sample,i, ido){
      df3 <- df_full[sample,]
      test_set <- !(c(1:nrow(df_full)) %in% sample)
     # print(ido)
      set.seed(42)
      rf_model_cv <- randomForest::randomForest(y=df3[,ido], x=df3[,variables],data = df3, importance=T,
                                                proximity =T, mtry=rfGrid[i,1], ntree = rfGrid[i,2])
      preds <- stats::predict(rf_model_cv, df_full[test_set,variables], type='response')
      rms <- Metrics::rmse(actual=df_full[test_set,ido], predicted=preds)
      co <- cor(df_full[test_set,ido],preds, method = 'pearson' )
      return(c(rms, co))
    }
    
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
    score_list <- lapply(samples_list, FUN = cv_rf, i=i, ido=id)
    
    rmse_list <- NULL
    cor_list <- NULL
    for (vector in score_list){
      rmse_list <- append(rmse_list,vector[1])
      cor_list <- append(cor_list,vector[2])
    }
    RMSEs <- mean(rmse_list, na.rm=T)
    CORs <- mean(cor_list, na.rm=T)
    if (flag == FALSE || is.null(rf_model)) {
      RMSEs <- NA
      CORs <- NA
    }
    optimisation_rf <- list(RMSEs, CORs, rf_model)
  }
  
  m_core <- NULL
  m_rmse <- NULL
  mods <- NULL
  for (i in 1:dim(rfGrid)[1]){ #
    set.seed(42)
   # print(id)
    vec <- optimisation_rf(i, rfGrid = rfGrid , variables=variables, df_full = df_full, id=id)
    m_rmse <- append(m_rmse, vec[[1]])
    m_core <- append(m_core, vec[[2]])
    mods[[i]] <- vec[[3]]
  }
  # Finding the parameters combination for which the mean auc of the cross-validation is maximized
  g <- which.min(m_rmse)
  if (length(g)>0){
    best_model <- c(rfGrid[g,1], rfGrid[g,2],  m_rmse[g], m_core[g])
  } else{
    best_model <- c(NA,NA,NA, NA)
  }
  # } else{
  #   best_model <- c(id, NA,NA,NA,NA, NA, NA, NaN)
  # }
  #write(id, paste('follow_nn_',,'.txt', sep=''), append=T)
  b_mod <- mods[[g]]
  return(list(best_model, b_mod))
}
  
