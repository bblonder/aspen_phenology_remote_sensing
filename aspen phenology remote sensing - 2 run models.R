library(terra)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ranger)
library(pdp)
library(viridis)
library(tidyr)
library(caret)

vars_phenology = c("OGI", "OGMn", "GSL.50", "EVImax")
df_all = read.csv('df_all.csv')


# make models
df_all_for_rf = df_all %>% 
  mutate(x_grid = cut(x, breaks=seq(min(x),max(x),by=100))) %>%
  mutate(y_grid = cut(y, breaks=seq(min(y),max(y),by=100))) %>%
  mutate(grid_id = as.numeric(factor(paste(x_grid, y_grid)))) %>%
  mutate(year=factor(as.character(year))) %>%
  mutate(cytotype=factor(cytotype,levels=c("diploid","triploid")))


split_train_test <- function(df, fraction, xvar, yvar)
{
  n_cells = max(df$grid_id)
  train_ids = sample(1:n_cells, fraction*n_cells, replace=FALSE)
  test_ids = setdiff(1:n_cells, train_ids)
  
  train=df %>% filter(grid_id %in% train_ids) %>% dplyr::select(any_of(c(xvar, yvar))) %>% na.omit
  test=df %>% filter(grid_id %in% test_ids) %>% dplyr::select(any_of(c(xvar, yvar))) %>% na.omit
  
  train_x = as.data.frame(train[,xvar])
  names(train_x) = xvar
  train_y = as.numeric(train[,yvar,drop=TRUE])
  test_x = as.data.frame(test[,xvar])
  names(test_x) = xvar
  test_y = as.numeric(test[,yvar,drop=TRUE])
  
  return(list(train_x=train_x, train_y=train_y, test_x=test_x, test_y=test_y))
}

train_test_model <- function(tt, xvar, categorical)
{
  print('training model')
  if (categorical==TRUE)
  {
    m_rf = ranger(x=tt$train_x, y=tt$train_y,
                  max.depth = 10,
                  num.trees = 1000,
                  classification=TRUE,
                  verbose=TRUE,
                  probability = TRUE,
                  write.forest=TRUE,
                  importance = "impurity")
  }
  else
  {
    m_rf = ranger(x=tt$train_x, y=tt$train_y,
                  max.depth = 10,
                  num.trees = 1000,
                  classification=FALSE,
                  verbose=TRUE,
                  probability = FALSE,
                  write.forest=TRUE,
                  importance = "impurity")    
  }
  
  print('testing model')
  labels_predicted = predict(m_rf, tt$test_x)$predictions
  
  print('summarizing model')
  if (categorical==TRUE)
  {
    # apply threshold to probabilistic model
    labels_predicted = as.numeric(factor(labels_predicted[,1] < 0.5))
    labels_true = as.numeric(factor(tt$test_y))
    
    cm = confusionMatrix(factor(labels_predicted), factor(labels_true))$byClass
    result_this = data.frame(t(cm), xvars=paste(xvar,collapse="*"), np_train = nrow(tt$train_x),nv=ncol(tt$train_x))
  }
  else
  {
    labels_true = tt$test_y
    
    stats = postResample(labels_true, labels_predicted)
    
    result_this = data.frame(rmse=stats["RMSE"], r2=stats["Rsquared"], xvars=paste(xvar,collapse="*"), np_train = nrow(tt$train_x),nv=ncol(tt$train_x))
  }
  
  return(list(result=result_this,model=m_rf))
}


try_models <- function(df, yvar, xvar, categorical=TRUE, iter=5, fraction=0.8)
{
  print(yvar)
  result = NULL
  models = vector(mode='list',length=iter)
  datasets = vector(mode='list',length=iter)
  for (i in 1:iter)
  {
    start_time <- Sys.time()
    
    print('splitting data')
    
    tt = split_train_test(df=df, xvar=xvar, yvar=yvar, fraction=fraction)
    
    output = train_test_model(tt, xvar, categorical=categorical)
    
    # copy over outputs
    result_this = output$result
    models[[i]] = output$model
    datasets[[i]] = tt
    
    # report last stats
    end_time <- Sys.time()
    result_this$runtime = end_time - start_time
    result_this$yvar = yvar
    
    result = rbind(result, result_this)
    
    print(i/iter)
  }
  
  return(list(result=result,models=models,datasets=datasets,xvar=xvar))
}

models_all = lapply(vars_phenology, try_models, 
                    df=df_all_for_rf,
                    categorical=FALSE,
                    iter=3,
                    fraction=0.8,
                    xvar=setdiff(names(df_all), c('aspen_cover','x','y',vars_phenology)))
names(models_all) = vars_phenology






















do_pdp <- function(model_list, pred.vars, grid.resolution, categorical, df_train)
{
  result = lapply(1:length(model_list$models), function(id) {
    print(id/length(model_list$models))
    
    df_train_this = na.omit(df_train[,model_list$xvar])
    
    result_this = partial(model_list$models[[id]], pred.var=pred.vars, 
                          grid.resolution=grid.resolution,
                          train=df_train_this,
                          prob=categorical,
                          progress='text')
    result_this$replicate=id
    
    return(result_this)
    
  })  
  
  result_all = do.call("rbind",result)
  return(result_all)
}

summarize_pdp <- function(pdp)
{
  yvar = pdp$yvar[1]
  
  pdp = pdp %>% as_tibble %>% dplyr::select(-yvar)
  
  if ("year" %in% names(pdp))
  {
    pdp = pdp %>% mutate(year = as.numeric(as.character(year))) 
  }
  if ("cytotype" %in% names(pdp))
  {
    pdp = pdp %>% mutate(cytotype = as.numeric(cytotype))
  }
  
  # need to do approximate groupings as the partial picks different predictor values based on samples...
  pdp_list = pdp %>% 
    group_by(replicate) %>% 
    group_split
  
  # stack up all the predictions by replicate into an array
  pdp_list = lapply(pdp_list, as.matrix)
  pdp_array = array(unlist(pdp_list), dim = c(nrow(pdp_list[[1]]), ncol(pdp_list[[1]]), length(pdp_list)))
  
  # calculate elementwise means
  q50 = apply(pdp_array, c(1,2), quantile, 0.50) %>% 
    as.data.frame
  names(q50) = names(pdp)
  
  q05 = apply(pdp_array, c(1,2), quantile, 0.05) %>% 
    as.data.frame
  names(q05) = names(pdp)
  q95 = apply(pdp_array, c(1,2), quantile, 0.95) %>% 
    as.data.frame
  names(q95) = names(pdp)
  
  # lose the replicate column
  final = q50 %>% 
    dplyr::select(-replicate, -yhat) %>%
    mutate(yhat.q50 = q50$yhat) %>%
    mutate(yhat.q05 = q05$yhat) %>%
    mutate(yhat.q95 = q95$yhat) %>%
    mutate(yvar=yvar)
  
  if ("cytotype" %in% names(pdp))
  {
    final = final %>% 
      mutate(cytotype = factor(cytotype,levels=c(1,2),labels=c("diploid","triploid"),ordered=TRUE))
  }
  
  return(final)
}


# make 1d pdps for each phenology response variable
pdps_all = lapply(models_all, make_all_pdps_1d)





save.image(file = 'workspace all.Rdata')




