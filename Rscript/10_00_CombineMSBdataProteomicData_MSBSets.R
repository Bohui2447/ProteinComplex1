# Combined MSB and proteomic data using MSB training/test sets

options(max.print = 100)
source("00_00_function_complex.R")
library(funcTools)

### MSB paper data, including training and testing data, with gene names
if (file.exists("../Rdata/10_00_trainingTestSets_featuresClean2.Rdata")){
  cat("load ../Rdata/10_00_trainingTestSets_featuresClean2.Rdata\n")
  load("../Rdata/10_00_trainingTestSets_featuresClean2.Rdata")
} else {
  load("../rawData/APMSdata/13trainingTestSets_features.Rdata")
  # All gene entrez IDs
  geneIDs = unique(c(trainingTestSets_featuresClean$geneid1, 
                     trainingTestSets_featuresClean$geneid2))
  
  # Make a data frame from IDs to gene names
  entrez2genename0 = ENTREZtoGeneName(geneIDs, sameOrder = F)
  entrez2genename1 = entrez2genename0[!duplicated(entrez2genename0$entrezgene),]
  entrez2genename = rbind(data.frame(entrezgene = setdiff(geneIDs, entrez2genename1$entrezgene),
                                     external_gene_name = "", stringsAsFactors = F), entrez2genename1)
  
  # Map IDs to gene names
  entrez2genename_clean1 = selectRows(entrez2genename, trainingTestSets_featuresClean$geneid1, 
                                      by = "entrezgene")$res
  entrez2genename_clean2 = selectRows(entrez2genename, trainingTestSets_featuresClean$geneid2, 
                                      by = "entrezgene")$res
  trainingTestSets_featuresClean$geneName1 = entrez2genename_clean1$external_gene_name
  trainingTestSets_featuresClean$geneName2 = entrez2genename_clean2$external_gene_name
  
  # remove unmapped rows
  id_notFound = (trainingTestSets_featuresClean$geneName1 == "" | trainingTestSets_featuresClean$geneName2 == "")
  trainingTestSets_featuresClean2 = trainingTestSets_featuresClean[!id_notFound,]
  View(head(trainingTestSets_featuresClean, 100))
  View(head(trainingTestSets_featuresClean2, 100))
  
  V14 = as.factor(c(trainP = 1, trainN = 0, testP = 1, testN = 0))
  trainingTestSets_featuresClean2$V14 = V14[as.character(trainingTestSets_featuresClean2$set)]
  
  # save
  rm(trainingTestSets_featuresClean, trainingTestSets_features, entrez2genename0, 
     entrez2genename1, entrez2genename_clean1, entrez2genename_clean2)
  save(trainingTestSets_featuresClean2, entrez2genename, 
       file = "../Rdata/10_00_trainingTestSets_featuresClean2.Rdata")
}

### According to gene names to obtain expression data
if (file.exists()){
  cat("load ../Rdata/10_00_exprAbs.Rdata\n")
  load("../Rdata/10_00_exprAbs.Rdata")
} else {
  idNames = trainingTestSets_featuresClean2[, c("geneid1", "geneid2", "geneName1", "geneName2", "V14")]
  
  # load PXD_combine dataset and complexes
  PXD_combine = get(load("../Rdata_localPC/3PXD_combine.Rdata"))
  # log data
  logExpr = log10(PXD_combine + 1)
  rm(PXD_combine, trainingTestSets_featuresClean2, entrez2genename); gc(T,T)
  # Absolute log changes
  exprAbs = abs(logExpr[idNames$geneName1,] - logExpr[idNames$geneName2,]) ## change to eucledian distance, to be done ***
  rownames(exprAbs) = rownames(idNames)
  rowNa = is.na(exprAbs[,1])
  
  save(idNames, exprAbs, rowNa, file = "../Rdata/10_00_exprAbs.Rdata")
}

# MSB features and absolute log changes combined
if(file.exists("../Rdata/10_00_featureExprAbs_pro1Pro2.Rdata")){
  cat("load ../Rdata/10_00_featureExprAbs_pro1Pro2.Rdata\n")
  load("../Rdata/10_00_featureExprAbs_pro1Pro2.Rdata")
} else {
  load("../Rdata/10_00_trainingTestSets_featuresClean2.Rdata")
  featureExprAbs_pro1Pro2 = cbind(trainingTestSets_featuresClean2, exprAbs)[!rowNa,]
  save(featureExprAbs_pro1Pro2, 
       file = "../Rdata/10_00_featureExprAbs_pro1Pro2.Rdata")
}

# Training, test, and validation sets  (preoteomic data)
if (file.exists("../Rdata/10_00_exprAllDivede_MSBTrain.Rdata") & 
    file.exists("../Rdata/10_00_exprAllDivede_MSBTest.Rdata") & 
    file.exists("../Rdata/10_00_exprAllDivede_MSBValidation.Rdata") & 
    file.exists("../Rdata/10_00_TrainMeanSD_MSBTrainTest.Rdata")){
  message("You have all training and testing and validation data, go to build models\n")
} else {
  # only absolute log changes (withou missing values)
  load("../Rdata/10_00_featureExprAbs_pro1Pro2.Rdata")
  
  # train, test and validation row index
  train_ind = which(featureExprAbs_pro1Pro2$set %in% c("trainP", "trainN"))
  testVal_ind = which(featureExprAbs_pro1Pro2$set %in% c("testP", "testN"))
  {
    set.seed(123)
    test_ind = sample(testVal_ind, size = floor(length(testVal_ind)/2), replace = F)
    val_ind = setdiff(testVal_ind, test_ind)
  }
  
  # output column index
  labelId = which(colnames(featureExprAbs_pro1Pro2) == "V14")
  
  ### training set
  trainset = as.matrix(featureExprAbs_pro1Pro2[train_ind, -(1:labelId)])
  { # mean and sd
    sdTrain = apply(trainset, 2, sd)
    meanTrain = colMeans(trainset)
  }
  
  noVarianceID = (sdTrain != 0)
  # normalization
  trainNorm = array(scale(trainset, center = meanTrain, scale = sdTrain), dim = dim(trainset))[,noVarianceID]
  rm(trainset); gc(T,T)# get more space
  # labels
  trainlables = as.factor(unname(featureExprAbs_pro1Pro2[train_ind, labelId]))
  trainlables = array(trainlables, dim = length(trainlables))
  # save data
  save(trainNorm, trainlables, 
       file =  "../Rdata/10_00_exprAllDivede_MSBTrain.Rdata")
  rm(trainNorm, trainlables); gc(T,T)
  
  ### test set
  testset = as.matrix(featureExprAbs_pro1Pro2[test_ind, -(1:labelId)])
  testNorm = array(scale(testset, center = meanTrain, scale = sdTrain), dim = dim(testset))[,noVarianceID]
  rm(testset); gc(T,T)# get more space
  testlables = as.factor(unname(featureExprAbs_pro1Pro2[test_ind, labelId]))
  testlables = array(testlables, dim = length(testlables))
  # save data
  save(testNorm, testlables, 
       file =  "../Rdata/10_00_exprAllDivede_MSBTest.Rdata")
  rm(testNorm, testlables); gc(T,T)
  
  ### validation set
  valset = as.matrix(featureExprAbs_pro1Pro2[val_ind, -(1:labelId)])
  valNorm = array(scale(valset, center = meanTrain, scale = sdTrain), dim = dim(valset))[,noVarianceID]
  rm(valset); gc(T,T)# get more space
  vallabels = as.factor(unname(featureExprAbs_pro1Pro2[val_ind, labelId]))
  vallabels = array(vallabels, dim = length(vallabels))
  # save data
  save(valNorm, vallabels, 
       file =  "../Rdata/10_00_exprAllDivede_MSBValidation.Rdata")
  rm(valNorm, vallabels); gc(T,T)
  
  rm(featureExprAbs_pro1Pro2); gc(T,T)
  
  # mean and variance
  meanTrain = meanTrain[noVarianceID]
  sdTrain = sdTrain[noVarianceID]
  save(meanTrain, sdTrain,
       file =  "../Rdata/10_00_TrainMeanSD_MSBTrainTest.Rdata")
  gc(T,T)
}




# Training, test, and validation sets (MSBdata + preoteomic data)
if (file.exists("../Rdata/10_00_MSBdataExpr_MSBTrain.Rdata") & 
    file.exists("../Rdata/10_00_MSBdataExpr_MSBTest.Rdata") & 
    file.exists("../Rdata/10_00_MSBdataExpr_MSBValidation.Rdata") & 
    file.exists("../Rdata/10_00_TrainMeanSD_MSBdataExpr_MSBTrainTest.Rdata")){
  message("You have all training and testing and validation data, go to build models\n")
} else {
  # only absolute log changes (withou missing values)
  load("../Rdata/10_00_featureExprAbs_pro1Pro2.Rdata")
  
  # train, test and validation row index
  train_ind = which(featureExprAbs_pro1Pro2$set %in% c("trainP", "trainN"))
  testVal_ind = which(featureExprAbs_pro1Pro2$set %in% c("testP", "testN"))
  {
    set.seed(123)
    test_ind = sample(testVal_ind, size = floor(length(testVal_ind)/2), replace = F)
    val_ind = setdiff(testVal_ind, test_ind)
  }
  
  # output column index
  labelId = which(colnames(featureExprAbs_pro1Pro2) == "V14")
  excludeId = which(colnames(featureExprAbs_pro1Pro2) %in% 
                      c("id", "set", "existHumap", "geneid1", "geneid2", 
                        "geneName1", "geneName2", "id1", "id2", "V14"))
  
  ### training set
  trainset = as.matrix(featureExprAbs_pro1Pro2[train_ind, -(excludeId)])
  if(FALSE){
  trainProPair_In1202 = rownames(trainset)
  save(trainProPair_In1202, file = '../Rdata/14_00_trainProPair_In1202.Rdata')
  }
  
  { # mean and sd
    sdTrain = apply(trainset, 2, sd)
    meanTrain = colMeans(trainset)
  }
  
  noVarianceID = (sdTrain != 0)
  # normalization
  trainNorm = array(scale(trainset, center = meanTrain, scale = sdTrain), dim = dim(trainset))[,noVarianceID]
  rm(trainset); gc(T,T)# get more space
  # labels
  trainlables = as.factor(unname(featureExprAbs_pro1Pro2[train_ind, labelId]))
  trainlables = array(trainlables, dim = length(trainlables))
  # save data
  save(trainNorm, trainlables, 
       file =  "../Rdata/10_00_MSBdataExpr_MSBTrain.Rdata")
  rm(trainNorm, trainlables); gc(T,T)# get more space
  
  ### test set
  testset = as.matrix(featureExprAbs_pro1Pro2[test_ind, -(excludeId)])
  if(FALSE){
    testProPair_In1202 = rownames(testset)
    save(testProPair_In1202, file = '../Rdata/14_00_testProPair_In1202.Rdata')
  }
  testNorm = array(scale(testset, center = meanTrain, scale = sdTrain), dim = dim(testset))[,noVarianceID]
  rm(testset); gc(T,T)# get more space
  # labels
  testlables = as.factor(unname(featureExprAbs_pro1Pro2[test_ind, labelId]))
  testlables = array(testlables, dim = length(testlables))
  # save data
  save(testNorm, testlables, 
       file =  "../Rdata/10_00_MSBdataExpr_MSBTest.Rdata")
  rm(testNorm, testlables); gc(T,T)# get more space
  
  ### validation set
  valset = as.matrix(featureExprAbs_pro1Pro2[val_ind, -(excludeId)])
  
  if(FALSE){
    valProPair_In1202 = rownames(valset)
    save(valProPair_In1202, file = '../Rdata/14_00_valProPair_In1202.Rdata')
  }
  
  valNorm = array(scale(valset, center = meanTrain, scale = sdTrain), dim = dim(valset))[,noVarianceID]
  rm(valset); gc(T,T)# get more space
  # labels
  vallabels = as.factor(unname(featureExprAbs_pro1Pro2[val_ind, labelId]))
  vallabels = array(vallabels, dim = length(vallabels))
  # save data
  save(valNorm, vallabels, 
       file =  "../Rdata/10_00_MSBdataExpr_MSBValidation.Rdata")
  rm(valNorm, vallabels); gc(T,T)# get more space
  
  rm(featureExprAbs_pro1Pro2); gc(T,T)
  # mean and variance
  meanTrain = meanTrain[noVarianceID]
  sdTrain = sdTrain[noVarianceID]
  save(meanTrain, sdTrain,
       file =  "../Rdata/10_00_TrainMeanSD_MSBdataExpr_MSBTrainTest.Rdata")
}
