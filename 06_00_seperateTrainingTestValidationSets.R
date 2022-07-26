# 'Rdata_localPC' is 'RdataAll'

rm(list = ls())

#### prepare data
if(file.exists("../Rdata_localPC/6exprAllDivede.Rdata") &
   file.exists("../Rdata_localPC/7TrainMeanSD.Rdata")){
  cat("load ../Rdata_localPC/6exprAllDivede.Rdata")
  load("../Rdata_localPC/6exprAllDivede.Rdata")
  
} else {
  # load("../Rdata_localPC/5exprLabel.Rdata")
  load("../Rdata_localPC/5exprLabel.Rdata")
  # exprLabel_norm = scale(exprLabel[, -ncol(exprLabel)])
  # exprLabel = as.matrix(cbind(exprLabel_norm, label = exprLabel$label))
  
  exprLabel = as.matrix(exprLabel)
  
  ## divide into train set and test set
  train_size = floor(0.95*nrow(exprLabel))
  # val_size = test_size = floor(0.05 * nrow(exprLabel))
  val_size = floor(0.025 * nrow(exprLabel))
  
  ## set index
  set.seed(123)
  train_ind = sort(sample(seq_len(nrow(exprLabel)), size = train_size, replace = FALSE))
  notInTrain = setdiff(seq_len(nrow(exprLabel)), train_ind)
  val_ind = sort(sample(notInTrain, size = val_size))
  test_ind = setdiff(notInTrain, val_ind)
  
  ## divide dataset
  trainset = exprLabel[train_ind, -ncol(exprLabel)]
  trainlables = as.factor(unname(exprLabel[train_ind, ncol(exprLabel)]))
  
  testset = exprLabel[test_ind, -ncol(exprLabel)]
  testlables = as.factor(unname(exprLabel[test_ind, ncol(exprLabel)]))
  
  valset = exprLabel[val_ind, -ncol(exprLabel)]
  vallabels = as.factor(unname(exprLabel[val_ind, ncol(exprLabel)]))
  
  #
  rm(exprLabel) # get more space for r
  
  ## normalization
  meanTrain = colMeans(trainset)
  sdTrain = apply(trainset, 2, sd)
  trainNorm = array(scale(trainset, center = meanTrain, scale = sdTrain), dim = dim(trainset))
  #
  rm(trainset) # get more space
  testNorm = array(scale(testset, center = meanTrain, scale = sdTrain), dim = dim(testset))
  valNorm = array(scale(valset, center = meanTrain, scale = sdTrain), dim = dim(valset))
  
  trainlables = array(trainlables, dim = length(trainlables))
  testlables = array(testlables, dim = length(testlables))
  vallabels = array(vallabels, dim = length(vallabels))
  
  # save data
  save(trainNorm, trainlables, testNorm, testlables, valNorm, vallabels, 
       file =  "../Rdata_localPC/6exprAllDivede.Rdata")
  save(meanTrain, sdTrain,
       file =  "../Rdata_localPC/7TrainMeanSD.Rdata")
}


