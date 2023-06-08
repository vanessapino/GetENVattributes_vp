PrepareDataset <- function(DataDirectory,CovariatesTarget,SpectraTarget,LocalDirectory,umap_,
                           ReferenceName,TargetProperty,TargetLocal,lower_waves,upper_waves,
                           neighbours_,COVARS,SpatialExtent){
  
  LocalCovars <- read.csv(paste0(DataDirectory,CovariatesTarget))
  
  
  LocalData <- read.csv(paste0(LocalDirectory,SpectraTarget))
  
  
  
  #Subset Global spectra based on desired wavelengths
  
  SpectralWaves <- colnames(GlobalDataset_Spectra)[grepl('scan_visnir',colnames(GlobalDataset_Spectra))]
  SpectralWaves_ <- sapply(SpectralWaves,function(x) gsub('scan_visnir.','',x))
  SpectralWaves_ <- as.numeric(sapply(SpectralWaves_,function(x) gsub('_pcnt','',x)))
  waves_to_use <-   SpectralWaves_>lower_waves&SpectralWaves_<upper_waves
  wavs <- SpectralWaves_[waves_to_use]
  
  #Check if  both global spectra and their respective observed values of Target Property have valid values
  
  IdswithValidSpec <- GlobalDataset_Spectra$id.layer_uuid_c[!is.na(GlobalDataset_Spectra$scan_visnir.2000_pcnt)]
  IdswithValidProp <- GlobalDataset_props$id.layer_uuid_c[!is.na(GlobalDataset_props[,TargetProperty])]
  
  SelectedSamples <- unique(IdswithValidProp,IdswithValidSpec)
  
  
  #Set Global spectral library with Target Property values and their respective environmental 
  # covariates in the exact same order 
  
  GlobalSpec <- GlobalSpec_All<- GlobalDataset_Spectra[match(GlobalDataset_Spectra$id.layer_uuid_c,SelectedSamples),SpectralWaves[waves_to_use]]
  GlobalCovars <- GlobalDataset[match(GlobalDataset$idlayer_uuid_c,SelectedSamples),]
  GlobalProp <- GlobalProp_All<- GlobalDataset_props[match(GlobalDataset_props$id.layer_uuid_c,SelectedSamples),TargetProperty]
  
  set.seed(124)
  toUmap <- GlobalCovars[complete.cases(GlobalCovars),COVARS]
  
  if(umap_){
    
    Projection <- uwot::tumap(toUmap,ret_model = T,n_neighbors = neighbours_)
    localProjected <- uwot::umap_transform(LocalCovars[,COVARS],Projection)
    
    rand_tr <-tri.mesh(localProjected[,1],localProjected[,2])
    rand.ch <- convex.hull(rand_tr,plot.it=F)
    pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))
    
    plot(Projection$embedding[,1],Projection$embedding[,2])
    
    points(localProjected[,1],localProjected[,2],col='Red',cex=5)
    lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue',lwd=5)
    inHull <- in.convex.hull(rand_tr,Projection$embedding[,1],Projection$embedding[,2])
    points(Projection$embedding[inHull,1],Projection$embedding[inHull,2],col='blue')
  }else{
    
    Projection <- prcomp(toUmap,center = T,scale. = T)
    localProjected <- predict(Projection,LocalCovars[,COVARS])
    
    rand_tr <-tri.mesh(localProjected[,1],localProjected[,2],duplicate = 'remove')
    rand.ch <- convex.hull(rand_tr,plot.it=F)
    pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))
    
    plot(Projection$x[,1],Projection$x[,2])
    points(localProjected[,1],localProjected[,2],col='Red')
    lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue')
    inHull <- in.convex.hull(rand_tr,Projection$x[,1],Projection$x[,2])
    points(Projection$x[inHull,1],Projection$x[inHull,2],col='blue')
  }
  
  #check in hull
  GlobalSpec_ <- GlobalSpec[inHull,]
  GlobalProp_ <- GlobalProp[inHull]
  GlobalCovars_ <- GlobalCovars[inHull,]
  
  # check spatial extent of library
  
  if (SpatialExtent){
    
    Locations_ <- lapply(GlobalCovars_$.geo[complete.cases(GlobalCovars_$.geo)],geojsonsf::geojson_sf)
    
    test <- Reduce(rbind,Locations_)
    
    sf::write_sf(test,paste0('Spatial/InHull_Locations_',ReferenceName,'_',neighbours_,'.shp'))
    
  }
  
  
  
  ##Create model with selected samples
  
  LocalSpectra <- LocalData[,grepl('^X[0-9]+',colnames(LocalData))]
  
  localwavs <- as.numeric(sapply(colnames(LocalSpectra),function(x) round(as.numeric(gsub('X','',x)))))
  colnames(LocalSpectra) <- sapply(colnames(LocalSpectra),function(x) paste0('scan_visnir.',round(as.numeric(gsub('X','',x))),'_pcnt'))
  
  #IF THE RESOLUTION IS NOT EXACTLY THE SAME AS THE GLOBAL DATASET#
  
  NewLocal <- prospectr::resample(LocalSpectra,wav = localwavs,new.wav = wavs)
  toPlotWavs <- as.numeric(colnames(NewLocal))
  colnames(NewLocal) <- sapply(colnames(NewLocal),function(x) paste0('scan_visnir.',round(as.numeric(gsub('X','',x))),'_pcnt'))
  
  
  ###Spectral Pre_Treatments
  
  GlobalSpec <- prospectr::detrend(GlobalSpec_,wav = wavs)
  GlobalSpec <- prospectr::savitzkyGolay(GlobalSpec,m = 0,p = 3,w = 11)
  GlobalSpec <- prospectr::standardNormalVariate(GlobalSpec)
  
  plot(toPlotWavs,NewLocal[1,],type='l',xlab = 'Wavelength nm',ylab = 'Reflectance')
  
  NewLocal <- prospectr::detrend(NewLocal,wav = wavs)
  
  plot(toPlotWavs,NewLocal[1,],type='l',xlab = 'Wavelength nm',ylab = 'Reflectance')
  
  NewLocal <- prospectr::savitzkyGolay(NewLocal,m = 0,p = 3,w = 11)
  plot(as.numeric(names(colnames(NewLocal))),NewLocal[1,],type='l',xlab = 'Wavelength nm',ylab = 'Reflectance')
  NewLocal <- prospectr::standardNormalVariate(NewLocal)
  plot(as.numeric(names(colnames(NewLocal))),NewLocal[1,],type='l',xlab = 'Wavelength nm',ylab = 'Reflectance')
  
  
  GlobalSpec_All <- prospectr::detrend(GlobalSpec_All,wav = wavs)
  GlobalSpec_All <- prospectr::savitzkyGolay(GlobalSpec_All,m = 0,p = 3,w = 11)
  GlobalSpec_All <- prospectr::standardNormalVariate(GlobalSpec_All)
  
  
  ###Variable Selection###
  
  ####
  
  ###Make model
  
  GlobalSpec <- GlobalSpec[complete.cases(GlobalProp_),]
  GlobalProp_ <- GlobalProp_[complete.cases(GlobalProp_)]
  print(paste0('model made with ',nrow(GlobalSpec),' samples'))
  
  GlobalSpec_All <- GlobalSpec_All[complete.cases(GlobalProp_All),]
  GlobalProp_All <- GlobalProp_All[complete.cases(GlobalProp_All)]
  
  
  ##Prepare for XGBoost
  
  Global_Training <- xgb.DMatrix(data = GlobalSpec_All,label=GlobalProp_All)
  Local_Training <- xgb.DMatrix(data = GlobalSpec,label=GlobalProp_)
  Test <- xgb.DMatrix(data = NewLocal)
  
  return(list(Global_Training=Global_Training,Local_Training=Local_Training,Test=Test))
}

RunXGB_ITER <- function(SEED=NULL,Global_Training,Local_Training,Test,neighbours_,Test_props){
  
  set.seed(SEED)
  params <- list(booster = "gbtree", 
                 objective = "reg:squarederror", 
                 eta=0.5, 
                 gamma=1, 
                 max_depth=6, 
                 min_child_weight=1, 
                 subsample=0.8, 
                 colsample_bytree=0.5)
  
  nrounds <- 50
  
  xgbcv_local <- xgb.cv( params = params,
                         data = Local_Training,
                         nrounds = nrounds,
                         nfold = 5,
                         showsd = T,
                         stratified = T,
                         print_every_n = 10,
                         early_stopping_rounds = 20,
                         maximize = F)
  
  
  xgb_Local <- xgb.train (params = params, 
                          data = Local_Training, 
                          nrounds = xgbcv_local$best_iteration, 
                          watchlist = list(train=Local_Training), 
                          print_every_n = 10, 
                          early_stopping_rounds = 10, 
                          maximize = F)
  
  xgbcv_Global <- xgb.cv( params = params,
                          data = Global_Training,
                          nrounds = nrounds,
                          nfold = 5,
                          showsd = T,
                          stratified = T,
                          print_every_n = 10,
                          early_stopping_rounds = 20,
                          maximize = F)
  
  xgb_Global <- xgb.train (params = params, 
                           data = Global_Training, 
                           nrounds = xgbcv_Global$best_iteration, 
                           watchlist = list(train=Global_Training), 
                           print_every_n = 10, 
                           early_stopping_rounds = 10, 
                           maximize = F )
  
  
  
  colnames(Test) <- NULL
  xgbpred_Local <- predict (xgb_Local,Test)
  
  
  
  Results <- rbind(data.frame(Type='Subset',Value=predict(xgb_Local,Test),X=Test_props$X,Y=Test_props$Y,UniqueID=Test_props$UniqueID),
                   data.frame(Type='Global',Value=predict(xgb_Global,Test),X=Test_props$X,Y=Test_props$Y,UniqueID=Test_props$UniqueID))
  
  print(plot(density(predict(xgb_Local,Test))))
  print(plot(density(predict(xgb_Global,Test))))
  
  
  xgb.save(xgb_Local,paste0('RData/Local_Model_',SEED,'_500_2300_','_n',neighbours_,'_',TargetLocal))
  xgb.save(xgb_Global,paste0('RData/Global_Model_',SEED,'_500_2300_','_n',neighbours_,'_',TargetLocal))
  
  print(ggplot(Results,aes(Value,col=Type,fill=Type))+
          geom_density(alpha=0.2)+
          theme_bw())
  
  return(Results)}
ReadXGB_ITER <- function(folder=NULL,property=NULL){
  
  models_ <- dir(folder,pattern = property,full.names = T)
  SEEDS <- unique(as.numeric(sapply(models_,function(x) stringr::str_extract(x,'(?<=Model_)[0-9]+'))))
  TestRMSEs <- do.call(rbind,lapply(SEEDS,function(model_){
    # model_ <- SEEDS[1]
    model_tmp <- models_[grepl(model_,models_)]
    
    Loc_model_ <- model_tmp[grepl('Local',model_tmp)]
    Glob_model_ <- model_tmp[grepl('Global',model_tmp)]
    xgb_Local <- xgb.load(Loc_model_)
    xgb_Global <- xgb.load(Glob_model_)
    
    xgbpred_Local <-  predict(xgb_Local,Test)
    xgbpred_Global <- predict(xgb_Global,Test)
    
    squared_error_Subset <- sqrt(mean((getinfo(Test,'label')-
                                         xgbpred_Local)^2))
    
    squared_error_Global <- sqrt(mean((getinfo(Test,'label')-
                                         xgbpred_Global)^2))
    # print(paste0('Iteration: ', model_,' Done'))
    return(rbind(data.frame(Type='Subset',Value=squared_error_Subset,Iter=model_),
                 data.frame(Type='Global',Value=squared_error_Global,Iter=model_)))
  }))
  
  return(TestRMSEs)
}  
ReadXGB_Preds <- function(folder=NULL,property=NULL){
  
  models_ <- dir(folder,pattern = property,full.names = T)
  SEEDS <- unique(as.numeric(sapply(models_,function(x) stringr::str_extract(x,'(?<=Model_)[0-9]+'))))
  Results <- do.call(rbind,lapply(SEEDS,function(model_){
    # model_ <- SEEDS[1]
    model_tmp <- models_[grepl(model_,models_)]
    
    Loc_model_ <- model_tmp[grepl('Local',model_tmp)]
    Glob_model_ <- model_tmp[grepl('Global',model_tmp)]
    xgb_Local <- xgb.load(Loc_model_)
    xgb_Global <- xgb.load(Glob_model_)
    
    xgbpred_Local <-  predict(xgb_Local,Test)
    xgbpred_Global <- predict(xgb_Global,Test)
    
    # print(paste0('Iteration: ', model_,' Done'))
    return(rbind(data.frame(Type='Subset',Value=xgbpred_Local,Iter=model_),
                 data.frame(Type='Global',Value=xgbpred_Global,Iter=model_)))
  }))
  
  return(Results)
}  
getRmse <- function(results_) {
  
  squared_error_Subset <- sqrt(mean((results_[results_$Type=='Observed','Value']-
                                       results_[results_$Type=='Subset','Value'])^2))
  
  squared_error_Global <- sqrt(mean((results_[results_$Type=='Observed','Value']-
                                       results_[results_$Type=='Global','Value'])^2))
  return(rbind(data.frame(Type='Subset',Value=squared_error_Subset),
               data.frame(Type='Global',Value=squared_error_Global)))
}
getShapImportance <- function(model_,trainset){
  shap_preps_local <- shap.prep(xgb_model = model_, X_train = trainset)
  
  waves_shap <- sapply(shap_preps_local$variable,function(x)as.numeric(gsub('_pcnt','',gsub('scan_visnir.','',x))))
  
  test <- lapply(unique(shap_preps_local$variable),function(x) {
    # browser()
    mean(shap_preps_local[grepl(x,shap_preps_local$variable),'value'][[1]],na.rm=TRUE)
  })
  
  tmp_ <- data.frame(Contribution=do.call(c,test),Variable=unique(waves_shap))
  ggplot(tmp_,aes(Variable,Contribution))+
    geom_line()
}
getShapValues <- function(model_,trainset){
  # model_ <- xgb_Local
  # trainset <-  NewLocal
  
  shap_preps_local <- shap.prep(xgb_model = model_, X_train = trainset)
  waves_shap <- as.character(unique(shap_preps_local$variable))
  
  test <- do.call(rbind,lapply(unique(shap_preps_local$ID),function(x) {
    # browser()
    tmp_ <- shap_preps_local[shap_preps_local$ID==x,]
    firstThree <- tmp_[order(tmp_$value,decreasing = T),'variable'][1:3][[1]]
    firstThree <- as.numeric(gsub('_pcnt','',gsub('scan_visnir.','',firstThree)))
    data.frame(ID=x,First=firstThree[1],
               Second=firstThree[2],
               Thirf=firstThree[3])
  }))
  
}


#   