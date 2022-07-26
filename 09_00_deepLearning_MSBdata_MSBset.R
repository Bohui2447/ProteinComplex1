# Using MSB training and test data to train a Deep learning model
# 
if (FALSE){
  
  library(funcTools)
  pathFile = rstudioapi::getActiveDocumentContext()
  path = dirname(pathFile$path)
  Rfile = pathFile$path
  
  email = "w.tao-2@umcutrecht.nl"
  set.seed(1234)
  param1 = seed = sample(1:1000000, 20, replace = F)
  logFolder = "09_00_logFileFolder"
  if(!dir.exists(paste0(path, "/", logFolder))){
    dir.create(paste0(path, "/", logFolder))
  }
  
  logFile = paste0(path, "/", logFolder, "/09_00_logfile_seed_", seed, ".log")
  res = rQsub(path = path, rFile = Rfile, jobName = paste0("09_00_jobSeed_", seed),
              threaded = 1, rTimeHour = 48, memoryG = 50,
              logFile = logFile, email = email, computeNode = "all.q@n00[3-7]*",
              preCMD = "echo \"module load R/3.5.1 && Rscript ",
              param1 = param1)
  res
  View(qstat())
  
}
# ### arguments 
args <- commandArgs(TRUE)
seed = as.numeric(args[1]) # random seed

message(seed)
.libPaths("/home/uu_npc/bli2/R/x86_64-pc-linux-gnu-library/3.5/")

message("Host")
message(system("echo $SGE_O_HOST"))

load("../Rdata/0_training_testing_features.Rdata")
source("00_00_function_complex.R")

# training data 
x_train = array(as.matrix(training[,-ncol(training)]), dim = dim(training[,-ncol(training)]))
y_train = array(c(N = '0', P = '1')[as.character(training[,ncol(training)])], dim = nrow(training))

# test data
x_test = array(as.matrix(testing[,-ncol(testing)]), dim = dim(testing[,-ncol(testing)]))
y_test = array(c(N = '0', P = '1')[as.character(testing[,ncol(testing)])], dim = nrow(testing))

library(keras)
library(funcTools)
library(ggplot2)
library(ROCR) # for roc.plot

print(seed)
set.seed(seed)
N = 100


params = list(dropOut1 = round(runif(N, 0, 0.5), digits = 3),
              dropOut2 = round(runif(N, 0, 0.5), digits = 3),
              dropOut3 = round(runif(N, 0, 0.5), digits = 3),
              units1 = sample(50*(4:15), size = N, replace = T),
              units2 = sample(10*(5:20), size = N, replace = T),
              units3 = sample(5*(2:10), size = N, replace = T))
for(ni in 1:N){
  message(ni)
  print(ni)
  {
    model = keras_model_sequential()  %>%
      layer_dense(units = params$units1[ni], input_shape = dim(x_train)[2], activation = "relu") %>% # 512
      layer_dropout(params$dropOut1[ni]) %>%
      layer_dense(units = params$units2[ni], activation = "relu") %>%
      layer_dropout(params$dropOut2[ni]) %>%
      layer_dense(units = params$units3[ni], activation = "relu") %>% 
      layer_dropout(params$dropOut3[ni]) %>%
      layer_dense(units = 1, activation = "sigmoid") %>% 
      compile(
        optimizer = optimizer_rmsprop(lr = 0.002),
        loss = 'binary_crossentropy',
        metrics = c('accuracy')
      )
    summary(model)
  }
  
  ## train the model
  set.seed(123)
  history = model %>% fit(x_train, as.numeric(y_train), #epochs=30, batch_size=10000, 
                          epochs = 50,
                          verbose = 1,
                          validation_data = list(x_test, as.numeric(y_test)))
  
  ## prediction for test set
  pre_prob = model %>% predict(x_test) %>% as.vector()
  pred <- prediction(predictions = pre_prob, labels = as.factor(y_test))
  
  ##### precision and recall
  prec_rec = performance(pred, "prec", "rec")
  prec_rec_plot = plot.prc(prec_rec)
  
  if (!dir.exists("../figs/09_00/precisionRecall")){
    dir.create("../figs/09_00/precisionRecall", recursive = T)
  }
  fileName = paste0("../figs/09_00/precisionRecall/09_00_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  message(fileName)
  ggsave(paste0(fileName, "_precisionRecall.png"),
         plot = prec_rec_plot, 
         width = 5, height = 5,dpi = 300)
  
  ##### ROC
  roc <- performance(pred,"tpr","tnr")
  auc <- performance(pred,"auc")
  auc = round(auc@y.values[[1]], 2)
  if (!dir.exists("../figs/09_00/ROC")){
    dir.create("../figs/09_00/ROC", recursive = T)
  }
  fileName = paste0("../figs/09_00/ROC/09_00_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  message(fileName)
  plotSave(paste0(fileName, "_ROC.png"),
           plotCMD = {
             plot(roc, colorize=TRUE)
             abline(a = 1, b = -1)
             legend(0, 0.2, auc, title = "AUC")
           }, width = 5, height = 5, dpi = 300)
  bestPerformance = cbind(prec_rec_plot$data[which.max(prec_rec_plot$data$F1),], AUC = auc)
  rownames(bestPerformance) = paste0("09_00_model_", ni, 
                                     "_seed", seed,"_",
                                     params$units1[ni], "_", params$dropOut1[ni],"_", 
                                     params$units2[ni], "_", params$dropOut2[ni], "_", 
                                     params$units3[ni],"_", params$dropOut3[ni])
  if (!dir.exists("../Rdata/09_00_bestPerformance")){
    dir.create("../Rdata/09_00_bestPerformance")
  }
  save(bestPerformance, file = paste0("../Rdata/09_00_bestPerformance/09_00_model_", ni, 
                                      "_seed", seed,"_",
                                      params$units1[ni], "_", params$dropOut1[ni],"_", 
                                      params$units2[ni], "_", params$dropOut2[ni], "_", 
                                      params$units3[ni],"_", params$dropOut3[ni], ".Rdata"))
  
  ##### plot history
  if (!dir.exists("../figs/09_00/history")){
    dir.create("../figs/09_00/history", recursive = T)
  }
  fileName = paste0("../figs/09_00/history/09_00_model_", ni, 
                    "_seed", seed,"_",
                    params$units1[ni], "_", params$dropOut1[ni],"_", 
                    params$units2[ni], "_", params$dropOut2[ni], "_", 
                    params$units3[ni],"_", params$dropOut3[ni])
  message(fileName)
  pl = plot(history)
  ggsave(paste0(fileName, "_history.png"), plot = pl, width = 7, height = 5,
         units = 'in', dpi = 300)
}

# # augment
# id = y_train[,1] == 2
# x_train_2 = x_train[id,]
# y_train_2 = y_train[id,]
# 
# n_pos = nrow(y_train)/sum(id)
# x_train_aug = rbind(x_train[!id,], matrixRep(x_train_2, m = floor(n_pos)))
# y_train_aug = rbind(y_train[!id,,drop = F], matrixRep(y_train_2, m = floor(n_pos)))
# 
# history = model %>% fit(x_train_aug, y_train_aug, epochs=30, batch_size=10000, 
#                         # x_test = x_test, y_test = x_test)
#                         validation_split = 0.2)
# plot(history)
# 
# 
# 
# # down sampling
# id = y_train[,1] == 2
# x_train_2 = x_train[id,]
# y_train_2 = y_train[id,,drop = FALSE]
# 
# id_1 = sample(1:sum(!id), size = sum(id), replace = F)
# x_train_down = rbind(x_train[!id,][id_1,], x_train_2)
# y_train_down = rbind(y_train[!id,,drop = F][id_1,,drop = F], y_train_2)
# 
# history = model %>% fit(x_train_down, y_train_down, epochs=30, 
#                         # x_test = x_test, y_test = x_test)
#                         validation_split = 0.2)
# plot(history)
# 

# # evaluate
# pre_prob = model %>% predict(x_test) %>% as.vector()
# pred <- prediction(predictions = pre_prob, labels = as.factor(y_test))
# prec_rec = performance(pred, "prec", "rec")
# prec_rec_plot = plot.prc(prec_rec)
# print(prec_rec_plot)

