# 'Rdata_localPC' is '2Rdata'

library(keras)
library(dplyr)
library(ggplot2)
library(purrr)
library(ROCR) # for roc.plot
library(funcTools)
setwd2thisFile()

rm(list = ls())
gc(reset = T)
#### load data
load("../Rdata_localPC/6exprAllDivede_ratio.Rdata")
source("00_00_function_complex.R")

#### Build the model
set.seed(123)
dropOut1 = round(runif(200, 0, 0.5), digits = 3)
dropOut2 = round(runif(200, 0, 0.5), digits = 3)
dropOut3 = round(runif(200, 0, 0.5), digits = 3)
units1 = sample(100*(1:7), size = 200, replace = T)
units2 = sample(50*(1:7), size = 200, replace = T)
units3 = sample(10*(1:7), size = 200, replace = T)
Model_performance = list()
for(ni in 1:200){
  print(ni)
  {
    model = keras_model_sequential()  %>%  
      layer_dense(units = units1[ni],  input_shape = dim(testNorm)[2], activation = "relu") %>% # 512
      # layer_dropout(0.3) %>%
      layer_dropout(dropOut1[ni]) %>%
      layer_dense(units = units2[ni], activation = "relu") %>%
      # layer_dropout(0.2) %>%
      layer_dropout(dropOut2[ni]) %>%
      layer_dense(units = units3[ni], activation = "relu") %>% # 64
      # layer_dropout(0.1) %>%
      layer_dropout(dropOut3[ni]) %>%
      layer_dense(units = 1, activation = "sigmoid")
    
    
    model  %>%  summary(model)
    
    model %>% compile(
      # optimizer_rmsprop(lr = 0.002, rho = 0.9, epsilon = NULL, decay = 0,
      # clipnorm = NULL, clipvalue = NULL),
      # optimizer = "rmsprop",
      optimizer = optimizer_rmsprop(lr = 0.001, rho = 0.9, decay = 0.01),
      loss = "binary_crossentropy",
      metrics = list('accuracy')
    )
    
    ## train the model
    set.seed(123)
    history = model %>% fit(
      # trainNorm, trainlables,
      trainNorm, as.numeric(trainlables),
      epochs = 10,
      verbose = 1,
      batch_size = 2048,
      # validation_data = list(valNorm, vallabels)
      validation_data = list(valNorm, as.numeric(vallabels))
    )
  }  
  # save performance
  Model_performance[[ni]] = history
  # save model
  fileName = paste0("../Rdata/10_01_model_", ni, "_", 
                    units1[ni], "_", dropOut1[ni],"_", 
                    units2[ni], "_", dropOut2[ni], "_", 
                    units3[ni],"_", dropOut3[ni])
  print(fileName)
  modelFileName = paste0(fileName,".h5")
  model %>% save_model_hdf5(modelFileName)
  
  
  ###### ROC curve
  pre_prob = model %>% predict(valNorm) %>% as.vector()
  pred <- prediction(predictions = pre_prob, labels = as.numeric(vallabels))
  
  # precision and recall
  prec_rec = performance(pred, "prec", "rec")
  # plot(prec_rec)
  prec_rec_plot = plot.prc(prec_rec)
  plot(prec_rec_plot)
  
  fileName = paste0("../figs/10_01_model_", ni, "_", 
                    units1[ni], "_", dropOut1[ni],"_", 
                    units2[ni], "_", dropOut2[ni], "_", 
                    units3[ni],"_", dropOut3[ni])
  ggsave(paste0(fileName, "_precisionRecall.png"),
         plot = prec_rec_plot, 
         width = 5, height = 5, res = 300)
  prec_rec_plot$data[which.max(prec_rec_plot$data$F1),]
  
  # ROC
  # roc <- performance(pred,"tpr","fpr")
  roc <- performance(pred,"tpr","tnr")
  # AUC
  auc <- performance(pred,"auc")
  auc = round(auc@y.values[[1]], 2)
  cat("AUC is", auc)
  
  plotSave(paste0(fileName, "_ROC.png"),
           plotCMD = {
             plot(roc, colorize=TRUE)
             abline(a = 1, b = -1)
             legend(0, 0.2, auc, title = "AUC")
           }, width = 5, height = 5, res = 300)
  
  # plot history
  pl = plot(history)
  ggsave(paste0(fileName, "_history.png"), plot = pl, width = 7, height = 5,
         units = 'in', dpi = 300)
}
save(Model_performance, file =  '../Rdata/10_01_Model_performance.Rdata')

