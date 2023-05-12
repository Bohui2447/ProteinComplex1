# MSB sets and proteomic data, then build Deep learning model.

if(FALSE){
  library(funcTools)
  # num = length(param1)
  pathFile = rstudioapi::getActiveDocumentContext()
  path = dirname(pathFile$path)
  Rfile = pathFile$path
  
  email = "b.li1@uu.nl"
  set.seed(1234)
  param1 = seed = sample(1:1000000, 40, replace = F)
  # param1 = seed = sample(1:1000000, 1, replace = F)
  
  dir.create(paste0(path, "/", "11_02_logFileFolder"))
  logFile = paste0(path, "/", "11_02_logFileFolder", "/11_02_logfile_seed_", seed, ".log")
  
  # rTimeHour = 4
  # memoryG = 30
  # rTimeHour = 40, memoryG = 70,
  res = rQsub2(path = path, rFile = Rfile, jobName = paste0("jobSeed_", seed),
              threaded = 1, rTimeHour = 50, memoryG = 70,
              logFile = logFile, email = email, 
              preCMD = "module load R/3.5.1 && Rscript ",
              param1 = param1)
  res
  View(qstat2(users = 'bli2'))
}
message(getwd())
args <- commandArgs(TRUE)
seed = as.numeric(args[1])
message(seed)
.libPaths("/home/uu_npc/bli2/R/x86_64-pc-linux-gnu-library/3.5/")
message(.libPaths())
options(max.print = 100)
library(keras)
library(dplyr)
library(ggplot2)
library(purrr)
library(ROCR) # for roc.plot
library(funcTools)
source("00_00_function_complex.R")
load("../Rdata/10_00_exprAllDivede_MSBTrain.Rdata") 
# load("../Rdata/10_00_exprAllDivede_MSBTest.Rdata") 
load("../Rdata/10_00_exprAllDivede_MSBValidation.Rdata") 

#### Build the model

print(seed)
set.seed(seed)
N =  50
# N = 2 # 100
params = list(dropOut1 = round(runif(N, 0, 0.5), digits = 3),
              dropOut2 = round(runif(N, 0, 0.5), digits = 3),
              dropOut3 = round(runif(N, 0, 0.5), digits = 3),
              units1 = sample(50*(4:15), size = N, replace = T),
              units2 = sample(10*(5:20), size = N, replace = T),
              units3 = sample(5*(2:10), size = N, replace = T))
for(ni in 1:N){
  print(ni)
  {
    model = keras_model_sequential()  %>%
      layer_dense(units = params$units1[ni], input_shape = dim(trainNorm)[2], activation = "relu") %>% # 512
      layer_dropout(params$dropOut1[ni]) %>%
      layer_dense(units = params$units2[ni], activation = "relu") %>%
      layer_dropout(params$dropOut2[ni]) %>%
      layer_dense(units = params$units3[ni], activation = "relu") %>% 
      layer_dropout(params$dropOut3[ni]) %>%
      layer_dense(units = 1, activation = "sigmoid") %>% 
      compile(
        optimizer = optimizer_rmsprop(lr = 0.002),
        loss = "binary_crossentropy",
        metrics = list('accuracy')
      )
    model  %>%  summary(model)
    
    
    ## train the model
    set.seed(123)
    history = model %>% fit(
      trainNorm, as.numeric(trainlables),
      epochs = 10,
      # epochs = 5,
      verbose = 1,
      batch_size = 2048,
      validation_data = list(valNorm, as.numeric(vallabels))
    )
  }
  
  ###### ROC curve
  pre_prob = model %>% predict(valNorm) %>% as.vector()
  pred <- prediction(predictions = pre_prob, labels = as.factor(vallabels))
  
  # precision and recall
  prec_rec = performance(pred, "prec", "rec")
  prec_rec_plot = plot.prc(prec_rec)
  
  if (!dir.exists("../figs/11_02/precisionRecall")){
    dir.create("../figs/11_02/precisionRecall", recursive = T)
  }
  fileName = paste0("../figs/11_02/precisionRecall/11_02_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  ggsave(paste0(fileName, "_precisionRecall.png"),
         plot = prec_rec_plot, 
         width = 5, height = 5,dpi = 300)
  
  F1_recall_precision = prec_rec_plot$data[which.max(prec_rec_plot$data$F1),]
  F1_value = max(prec_rec_plot$data$F1)
  
  # ROC
  roc <- performance(pred,"tpr","tnr")
  # AUC
  auc <- performance(pred,"auc")
  auc = round(auc@y.values[[1]], 2)
  cat("AUC is", auc)
  
  if (!dir.exists("../figs/11_02/ROC")){
    dir.create("../figs/11_02/ROC", recursive = T)
  }
  fileName = paste0("../figs/11_02/ROC/11_02_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  plotSave(paste0(fileName, "_ROC.png"),
           plotCMD = {
             plot(roc, colorize=TRUE)
             abline(a = 1, b = -1)
             legend(0, 0.2, auc, title = "AUC")
           }, width = 5, height = 5, dpi = 300)
  
  # plot history
  if (!dir.exists("../figs/11_02/history")){
    dir.create("../figs/11_02/history", recursive = T)
  }
  fileName = paste0("../figs/11_02/history/11_02_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  pl = plot(history)
  ggsave(paste0(fileName, "_history.png"), plot = pl, width = 7, height = 5,
         units = 'in', dpi = 300)
  
  ## Save model and results
  if (!dir.exists("../Rdata/ModelSave/11_02/model/")){
    dir.create("../Rdata/ModelSave/11_02/model/", recursive = T)
  }
  # save model
  if (F1_value > 0.49){
  # if (F1_value > 0.4){
    fileName = paste0("../Rdata/ModelSave/11_02/model/11_02_model_", ni, 
                      "_seed", seed,"_",
                      params$units1[ni], "_", params$dropOut1[ni],"_", 
                      params$units2[ni], "_", params$dropOut2[ni], "_", 
                      params$units3[ni],"_", params$dropOut3[ni])
    model %>% save_model_hdf5(paste0(fileName,"_model.h5"))
    
    # save parameters
    if (!dir.exists("../Rdata/ModelSave/11_02/parameter/")){
      dir.create("../Rdata/ModelSave/11_02/parameter/", recursive = T)
    }
    fileName = paste0("../Rdata/ModelSave/11_02/parameter/11_02_model_", ni, 
                      "_seed", seed,"_",
                      params$units1[ni], "_", params$dropOut1[ni],"_", 
                      params$units2[ni], "_", params$dropOut2[ni], "_", 
                      params$units3[ni],"_", params$dropOut3[ni])
    pre_prob_trainSet = model %>% predict(trainNorm) %>% as.vector()
    pred_trainSet <- prediction(predictions = pre_prob_trainSet,
                                labels = as.factor(trainlables))
    
    parameter = list(F1_recall_precision = F1_recall_precision,
                     auc  = auc, roc = roc, pl = pl, 
                     prec_rec_plot = prec_rec_plot,
                     pre_prob = pre_prob, pred = pred, 
                     pre_prob_trainSet = pre_prob_trainSet,
                     pred_trainSet = pred_trainSet)
    save(parameter, file = paste0(fileName,'_Parameter.Rdata'))
  }
}
