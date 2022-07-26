# MSB sets, MSB data + proteomic data, then build Deep learning model. 
# 12_02 only use validation set to evalute model performance

## 12_03 plus test set ppis to re-caculate model performance





options(max.print = 100)
library(keras)
library(dplyr)
library(ggplot2)
library(purrr)
library(ROCR) # for roc.plot
library(funcTools)
source("00_00_function_complex.R")

## combine test and validation data sets
# load("../Rdata/10_00_MSBdataExpr_MSBTrain.Rdata") 
if(file.exists('../Rdata/12_03_MSBdataExpr_Combine_TextValidation.Rdata')){
  cat('load ../Rdata/12_03_MSBdataExpr_Combine_TextValidation.Rdata')
  load('../Rdata/12_03_MSBdataExpr_Combine_TextValidation.Rdata')
}else {
  
  load("../Rdata/10_00_MSBdataExpr_MSBTest.Rdata")
  load("../Rdata/10_00_MSBdataExpr_MSBValidation.Rdata") 
  # combine valnorm and test norm
  testNorm_combine = rbind(valNorm, testNorm)
  # combine label
  lableCombine = c(vallabels, testlables)
  test_combineLables = array(lableCombine, dim = length(lableCombine))
  # save data
  save(testNorm_combine, test_combineLables, 
       file = '../Rdata/12_03_MSBdataExpr_Combine_TextValidation.Rdata')
}



##### Recaculate model performance
load("../Rdata/10_00_MSBdataExpr_MSBTrain.Rdata") 
load('../Rdata/12_03_MSBdataExpr_Combine_TextValidation.Rdata')

folderBase = "12_03"

fl = list.files('../Rdata/ModelSave/12_02/model/')
for(ni in fl){
  print(ni)
  model = load_model_hdf5(paste0('../Rdata/ModelSave/12_02/model/', ni))
  
  ###### ROC curve
  pre_prob = model %>% predict(testNorm_combine) %>% as.vector()
  pred <- prediction(predictions = pre_prob, labels = as.factor(test_combineLables))
  
  # precision and recall
  prec_rec = performance(pred, "prec", "rec")
  prec_rec_plot = plot.prc(prec_rec)
  
  if (!dir.exists(paste0("../figs/",folderBase, "/precisionRecall"))){
    dir.create(paste0("../figs/", folderBase, "/precisionRecall"), recursive = T)
  }
  fileName = paste0("../figs/", folderBase, "/precisionRecall/" ,ni)
  ggsave(paste0(fileName, "_precisionRecall.png"),
         plot = prec_rec_plot, 
         width = 5, height = 5,dpi = 300)
  F1_recall_precision = prec_rec_plot$data[which.max(prec_rec_plot$data$F1),]
  
  F1_value = max(prec_rec_plot$data$F1)
  print(F1_value)
  
  # ROC
  roc <- performance(pred,"tpr","tnr")
  # AUC
  auc <- performance(pred,"auc")
  auc = round(auc@y.values[[1]], 2)
  cat("AUC is", auc)
  
  if (!dir.exists(paste0("../figs/", folderBase, "/ROC"))){
    dir.create(paste0("../figs/", folderBase, "/ROC"), recursive = T)
  }
  fileName = paste0("../figs/", folderBase, "/ROC/", ni)
  plotSave(paste0(fileName, "_ROC.png"),
           plotCMD = {
             plot(roc, colorize=TRUE)
             abline(a = 1, b = -1)
             legend(0, 0.2, auc, title = "AUC")
           }, width = 5, height = 5, dpi = 300)
  
  # # plot history
  # if (!dir.exists(paste0("../figs/", folderBase, "/history"))){
  #   dir.create(paste0("../figs/", folderBase, "/history"), recursive = T)
  # }
  # fileName = paste0("../figs/", folderBase, "/history/", ni)
  # pl = plot(history)
  # ggsave(paste0(fileName, "_history.png"), plot = pl, width = 7, height = 5,
  #        units = 'in', dpi = 300)
  # 
  ## Save model and results
  if (!dir.exists(paste0("../Rdata/ModelSave/", folderBase, "/model/"))){
    dir.create(paste0("../Rdata/ModelSave/", folderBase, "/model/"), recursive = T)
  }
  # save model
  if (F1_value > 0.60){
    # if (F1_value > 0.5){
    fileName = paste0("../Rdata/ModelSave/", folderBase, "/model/", ni)
    model %>% save_model_hdf5(paste0(fileName,"_model.h5"))
    
    # save parameters
    if (!dir.exists(paste0("../Rdata/ModelSave/", folderBase, "/parameter/"))){
      dir.create(paste0("../Rdata/ModelSave/", folderBase, "/parameter/"), recursive = T)
    }
    fileName = paste0("../Rdata/ModelSave/", folderBase, "/parameter/", ni)
    pre_prob_trainSet = model %>% predict(trainNorm) %>% as.vector()
    pred_trainSet <- prediction(predictions = pre_prob_trainSet,
                                labels = as.factor(trainlables))
    
    parameter = list(F1_recall_precision = F1_recall_precision,
                     auc  = auc, roc = roc,  
                     prec_rec_plot = prec_rec_plot,
                     pre_prob = pre_prob, pred = pred, 
                     pre_prob_trainSet = pre_prob_trainSet,
                     pred_trainSet = pred_trainSet)
    save(parameter, file = paste0(fileName,'_Parameter.Rdata'))
  }
}
