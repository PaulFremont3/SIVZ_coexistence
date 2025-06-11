best_models_lm<-function(id){
  optimisation_gam <- function(i, gamGrid, variables, nsp, df_full){
    flag <- TRUE
    gam_model <- NULL
    # Fitting the model on the whole dataset
    set.seed(42)
    possibleError <- tryCatch(
      gam_model <- lm(y ~ log10Vhost + RVirus+VirusType+HostType ,
                             data = df_full),
      error=function(e) flag <- FALSE
    )
    cv_gam <- function(sample,i, nsp, df_full){
      df3 <- df_full[sample,]
      test_set <- !(c(1:nrow(df_full)) %in% sample)
      set.seed(42)
      gam_model_cv <- lm(y ~ log10Vhost   + RVirus+VirusType+HostType ,
                                data = df3)
      preds <- stats::predict(gam_model_cv, df_full[test_set,variables], type='response')
      #print(length(preds))
      #print(preds)
      rms <- Metrics::rmse(actual=10^df_full[test_set,id], predicted=10^preds)
      rms_bis <- Metrics::rmse(actual=df_full[test_set,id], predicted=preds)
      co <- cor(df_full[test_set,id],preds, method = 'pearson' )
      return(c(rms, co, rms_bis))
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
    score_list <- lapply(samples_list, FUN = cv_gam, i=i, nsp=nsp, df_full=df_full)
    
    rmse_list <- NULL
    rmse_list_bis <- NULL
    cor_list <- NULL
    for (vector in score_list){
      rmse_list <- append(rmse_list,vector[1])
      cor_list <- append(cor_list,vector[2])
      rmse_list_bis <- append(rmse_list_bis,vector[3])
    }
    RMSEs <- mean(rmse_list, na.rm=T)
    RMSEs_bis <- mean(rmse_list_bis, na.rm=T)
    CORs <- mean(cor_list, na.rm=T)
    #print(cor_list)
    if (flag == FALSE || is.null(gam_model)) {
      RMSEs <- NA
      RMSEs_bis <- NA
      CORs <- NA
    }
    optimisation_gam <- list(RMSEs, CORs, RMSEs_bis, gam_model)
  }
  
  set.seed(42)
  df_full <- as.data.frame(df_full)
  
  m_core <- NULL
  m_rmse <- NULL
  m_rmse_bis <- NULL
  mods <- list()
  for (i in 1:dim(gamGrid)[1]){
    set.seed(42)
    vec <- optimisation_gam(i, gamGrid = gamGrid , variables=variables, nsp = gamGrid[i,1], df_full=df_full)
    m_rmse <- append(m_rmse, vec[[1]])
    m_rmse_bis <- append(m_rmse_bis, vec[[3]])
    m_core <- append(m_core, vec[[2]])
    mods[[i]] <- vec[[4]]
  }
  #print(m_rmse)
  #print(m_core)
  # Finding the parameters combination for which the mean auc of the cross-validation is maximized
  g <- which.min(m_rmse_bis)
  if (length(g)>0){
    best_model <- c(gamGrid[g,1],  m_rmse[g], m_rmse_bis[g], m_core[g])
  } else{
    best_model <- c(NA,NA,NA, NA)
  }
  b_mod <- mods[[g]]
  print(id)
  return(list(best_model,b_mod ))
  
}
