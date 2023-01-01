# MSB sets, MSB data + proteomic data, then build Deep learning model. ### To be finished ********

if(FALSE){
  library(funcTools)
  # num = length(param1)
  pathFile = rstudioapi::getActiveDocumentContext()
  path = dirname(pathFile$path)
  Rfile = pathFile$path
  
  email = "w.tao-2@umcutrecht.nl"
  set.seed(1234)
  param1 = seed = sample(1:1000000, 20, replace = F)
  dir.create(paste0(path, "/", "12_00_logFileFolder"))
  logFile = paste0(path, "/", "12_00_logFileFolder", "/12_00_logfile_seed_", seed, ".log")
  res = rQsub(path = path, rFile = Rfile, jobName = paste0("jobSeed_", seed),
              threaded = 1, rTimeHour = 48, memoryG = 90,
              logFile = logFile, email = email, computeNode = "all.q@n00[3-7]*",
              preCMD = "echo \"module load R/3.5.1 && Rscript ",
              param1 = param1)
  res
  View(qstat())
}

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
load("../Rdata/10_00_MSBdataExpr_MSBTrain.Rdata") 
# load("../Rdata/10_00_MSBdataExpr_MSBTest.Rdata") 
load("../Rdata/10_00_MSBdataExpr_MSBValidation.Rdata") 

#### Build the model

print(seed)
set.seed(seed)
N = 100
params = list(dropOut1 = round(runif(N, 0, 0.5), digits = 3),
              dropOut2 = round(runif(N, 0, 0.5), digits = 3),
              dropOut3 = round(runif(N, 0, 0.5), digits = 3),
              units1 = sample(50*(4:15), size = N, replace = T),
              units2 = sample(10*(5:20), size = N, replace = T),
              units3 = sample(5*(2:10), size = N, replace = T))
for(ni in N:1){
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
      # epochs = 50,
      epochs = 10,
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
  
  if (!dir.exists("../figs/12_00/precisionRecall")){
    dir.create("../figs/12_00/precisionRecall", recursive = T)
  }
  fileName = paste0("../figs/12_00/precisionRecall/12_00_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  ggsave(paste0(fileName, "_precisionRecall.png"),
         plot = prec_rec_plot, 
         width = 5, height = 5,dpi = 300)
  prec_rec_plot$data[which.max(prec_rec_plot$data$F1),]
  
  # ROC
  roc <- performance(pred,"tpr","tnr")
  # AUC
  auc <- performance(pred,"auc")
  auc = round(auc@y.values[[1]], 2)
  cat("AUC is", auc)
  
  if (!dir.exists("../figs/12_00/ROC")){
    dir.create("../figs/12_00/ROC", recursive = T)
  }
  fileName = paste0("../figs/12_00/ROC/12_00_model_", ni, 
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
  if (!dir.exists("../figs/12_00/history")){
    dir.create("../figs/12_00/history", recursive = T)
  }
  fileName = paste0("../figs/12_00/history/12_00_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  pl = plot(history)
  ggsave(paste0(fileName, "_history.png"), plot = pl, width = 7, height = 5,
         units = 'in', dpi = 300)
}
