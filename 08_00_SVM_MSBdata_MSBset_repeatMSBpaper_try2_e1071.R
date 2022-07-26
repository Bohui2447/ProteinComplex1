if (FALSE){
  setwd("/hpc/uu_npc/bli2/ProteinComplex/Rscript/")
  load("../rawData/APMSdata/13trainingTestSets_features.Rdata")
  
  trainingTestSets_featuresClean$V14 = as.factor(c(trainP = 1, trainN = 0,
                                                   testP = 1, testN = 0))[trainingTestSets_featuresClean$set]
  
  # intrain <- createDataPartition(y = trainingTestSets_featuresClean$V14, p= 0.1, list = FALSE)
  # intest <- createDataPartition(y = trainingTestSets_featuresClean$V14, p= 0.1, list = FALSE)
  # 
  # training = subset(trainingTestSets_featuresClean[intrain[,1],], set %in% c("trainP", "trainN"))[, -(1:9)]
  # testing = subset(trainingTestSets_featuresClean[intest[,1],], set %in% c("testP", "testN"))[, -(1:9)]
  
  training0 = subset(trainingTestSets_featuresClean, set %in% c("trainP", "trainN"))[, -(1:9)]
  testing0 = subset(trainingTestSets_featuresClean, set %in% c("testP", "testN"))[, -(1:9)]
  
  
  rm(trainingTestSets_featuresClean, trainingTestSets_features)
  gc(T,T)
  
  training = scale(training0[, -ncol(training0)])
  ms = attributes(training)
  training = as.data.frame(training)
  training$V14 = training0$V14
  training$V14 = factor(c("N", "P"), c("N", "P"))[as.integer(training$V14)]
  
  testing = as.data.frame(scale(testing0[, -ncol(testing0)], center = ms$`scaled:center`, 
                                scale = ms$`scaled:scale`))
  testing$V14 = testing0$V14
  testing$V14 = factor(c("N", "P"), c("N", "P"))[as.integer(testing$V14)]
  
  save(training, testing, file = "../Rdata/0_training_testing_features.Rdata")
}

if (FALSE){
  library(funcTools)
  options(max.print = 100)
  
  # num = length(param1)
  # path = "/hpc/uu_npc/bli2/ProteinComplex/"
  # Rfile = paste0(path, "Rscript/7_runClusterOneOnScore.R")
  
  Rfile = rstudioapi::getActiveDocumentContext()$path
  path = dirname(Rfile)
  # logFile = paste0(Rfile, "_logFile.log")
  
  # logFile = paste0(path, "Rscript/5_logfilename.log")
  # email = "b.li1@uu.nl"
  email = "w.tao@umcutrecht.nl"
  
  # grid_radial <- expand.grid(sigma = c(0.01, 0.02, 0.025, 0.03, 0.04,
  #                                      0.05, 0.06, 0.07,0.08, 0.09, 0.1, 
  #                                      0.25, 0.5, 0.75,0.9),
  #                            C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
  #                                  1, 1.5, 2, 5))#[c(87, 123),]
  
  # grid_radial <- expand.grid(sigma = 2^seq(3, -15, by = -2),
  #                            C = 2^seq(-5, 15, by = 2))[(c(8, 63)+2),]
  grid_radial <- expand.grid(sigma = 2^seq(3, -15, by = -2),
                             C = 2^seq(-5, 15, by = 2))[-c(8:10, 64:65),]
  
  
  param1 = as.character(grid_radial[,1])
  param2 = as.character(grid_radial[,2])
  
  tmp = apply(grid_radial, 1, function(x)paste0(x, collapse = "_"))
  logFiles = paste0(path, "5_logfilename_", tmp, "_e1071.log")
  jobNames = paste0("e1071_", tmp)
  
  # qsub Rjobs
  # res = rQsub(path = path, rFile = Rfile, jobName = jobNames[1:2],
  #             threaded = 1, rTimeHour = 200, memoryG = 20,
  #             logFile = logFiles[1:2], email = email,
  #             preCMD = "echo \"module load R/3.4.1 && Rscript ",
  #             param1 = param1[1:2], param2 = param2[1:2], param3 = param3[1:2])
  
  res = rQsub(path = path, rFile = Rfile, jobName = jobNames,
              threaded = 1, rTimeHour = 200, memoryG = 30,
              logFile = logFiles, email = email,
              preCMD = "echo \"module load R/3.5.1 && Rscript ",
              param1 = param1, param2 = param2)
  
  res
  qstat()
}

### arguments 
args <- commandArgs(TRUE)
ar1 = as.numeric(args[1]) # sigma
ar2 = as.numeric(args[2]) # C

library(methods)
library(e1071)
setwd("/hpc/uu_npc/bli2/ProteinComplex/Rscript/")
load("../Rdata/0_training_testing_features.Rdata")


# grid_radial0 <- expand.grid(sigma = ar1, C = ar2)
set.seed(3233)
svmFit = tune.svm(V14 ~ ., data = training, 
                  gamma = ar1, 
                  cost = ar2, 
                  probability = T, 
                  tunecontrol = tune.control(cross = 10))
message(paste("Best parameters:", svmFit$best.parameters))
message(paste("Best performance:", svmFit$best.performance))

warnings()
save(svmFit, file = paste0("../Rdata/e1071_svmFit_",ar1, "_", ar2, "_cv10.Rdata"))

if (T){
  # plot(svm_Radial_Grid)
  svm_Radial_Grid = svmFit$best.model
  test_pred_Radial_Grid = attributes(predict(svm_Radial_Grid, testing[,(colnames(testing) != "V14")], 
                                             probability = T))$probabilities
}
test_pred_Radial_Grid = as.data.frame(test_pred_Radial_Grid)
save(svmFit, test_pred_Radial_Grid, file = paste0("../Rdata/e1071_svmFit_",
                                                  ar1, "_", ar2, "_cv10.Rdata"))

if (FALSE){
  
  
  library(ggplot2)
  plot.prc = function(prc.perf, main = NULL, xlim = c(0, 1), ylim = c(0,1), 
                      plotF1 = TRUE, nudge_x = 0.02){
    F1 = function(x) 2*x[1]*x[2]/(x[1]+x[2])
    
    data_prc = data.frame(prc.perf@x.values, prc.perf@y.values)
    colnames(data_prc) = c(prc.perf@x.name, prc.perf@y.name)
    if(is.na(data_prc[1,1])){
      data_prc[1,1] = 0
    } 
    if(is.na(data_prc[1,2])){
      data_prc[1,2] = 1
    } 
    data_prc$F1 = apply(data_prc[,1:2], 1, F1)
    id = which.max(data_prc$F1)[1]
    data_F1 = data_prc[id,]
    g = ggplot(data_prc, aes_string(x = prc.perf@x.name, y = prc.perf@y.name))+
      geom_line()+
      geom_point(data = data_F1)+
      coord_cartesian(xlim = xlim, ylim = ylim)+
      ggtitle(main)+
      theme_Publication()
    if(plotF1) g = g +
      geom_text(aes(label = paste0("F1=", round(F1, 3))), data = data_F1, 
                hjust = "outward", nudge_x = nudge_x, color = "brown")
    return(g)
  }
  
  # https://www.r-bloggers.com/support-vector-machine-simplified-using-r/
  # test_pred_Radial_Grid <- predict(svm_Radial_Grid, newdata = testing, type = "prob")
  library(ROCR)
  load("../Rdata/0_training_testing_features.Rdata")
  resFiles = grep("e1071_svmFit_", dir("../Rdata/", full.names = T), value = T)
  for (mi in resFiles){
    print(mi)
    load(mi)
    print(head(test_pred_Radial_Grid))
    if(is.matrix(test_pred_Radial_Grid)){
      test_pred_Radial_Grid = as.data.frame(test_pred_Radial_Grid)  
    }
    if (!all(is.na(test_pred_Radial_Grid$P))){
      pred <- prediction(test_pred_Radial_Grid$P, testing$V14)
      # roc.perf = performance(pred, measure="prec", x.measure="rec")
      prc.perf = performance(pred, measure="prec", x.measure="rec")
      # plot(roc.perf, colorize = T)
      figFile = strSplit(mi, split = "/")[, 4]
      
      # plotSave(filename = paste0("../figs/", figFile, ".png"), 
      #          plotCMD = plot(prc.perf, main = figFile, xlim = c(0, 1), ylim = c(0,1)),
      #          width = 6.5, height = 5)
      plotSave(filename = paste0("../figs/", figFile, ".png"), 
               plotCMD = plot.prc(prc.perf, main = figFile, xlim = c(0, 1), ylim = c(0,1)),
               width = 6.5, height = 5)
    }
  }
}

