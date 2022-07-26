
library(funcTools)
library(keras)
source('00_00_function_complex.R')
options(max.print = 50)


######## use protein pairs from MSB to broden unknown (random) ppis, and then predict ppi score using the best model
# ppis with ap-ms features from MSB
if(file.exists('../Rdata/14_00_ppi_predictProb_all.Rdata')){
  cat('load ../Rdata/14_00_ppi_predictProb_all.Rdata')
  load('../Rdata/14_00_ppi_predictProb_all.Rdata')
} else {
  humapFeature = get(load('../Rdata/15_00_humapFeature_transformGeneidToGenesymble.Rdata'))
  humapFeature = humapFeature[, c(1:6, 265, 266, 7:264)]
  # load best model
  model = load_model_hdf5('../Rdata/ModelSave/12_03/model/12_02_model_18_seed186720_350_0.438_140_0.214_25_0.037_model.h5_model.h5')
  # load mean and sd of training dataset
  load('../Rdata/10_00_TrainMeanSD_MSBdataExpr_MSBTrainTest.Rdata')
  
  
  ### combine ap-ms feature to ms/ms feature
  PXD_combine = get(load("../Rdata_localPC/3PXD_combine.Rdata"))
  logExpr = log10(PXD_combine + 1)
  
  ## split humapFeature into small dataset to predict ppi
  
  ppi_predictProb = list()
  i1 = seq(1, 2518356, by = 100000)
  i2 = c(i1[-1]-1, 2518356)
  for(ni in 1:length(i1)){
    
    expid1 = logExpr[humapFeature$geneSymble1[i1[ni]:i2[ni]],]
    expid2 = logExpr[humapFeature$geneSymble2[i1[ni]:i2[ni]],]
    exprAbs = abs(expid1 - expid2)
    rm(expid1, expid2);gc(T,T)
    apmsExpr = cbind(humapFeature[i1[ni]:i2[ni], 9:ncol(humapFeature)], exprAbs)
    rm(exprAbs);gc(T,T)
    # remove columns that sd = 0
    apmsExpr = apmsExpr[, -(which(!colnames(apmsExpr) %in% names(sdTrain)))]
    # normalization
    apmsExprNorm = array(scale(apmsExpr, center = meanTrain, scale = sdTrain), dim = dim(apmsExpr))
    rm(apmsExpr); gc(T,T)
    ## predict ppi 
    time = system.time({
      pre_prob = model %>% predict(apmsExprNorm) %>% as.vector()
    })
    print(time)
    rm(apmsExprNorm); gc(T,T)
    # pred <- prediction(predictions = pre_prob, labels = as.factor(vallabels))
    pre_prob = data.frame(geneSymble1 = humapFeature$geneSymble1[i1[ni]:i2[ni]],
                          geneSymble2 = humapFeature$geneSymble2[i1[ni]:i2[ni]],
                          pre_prob)
    rownames(pre_prob) = paste0(humapFeature$geneid1[i1[ni]:i2[ni]], "_", 
                                humapFeature$geneid2[i1[ni]:i2[ni]])
    ppi_predictProb[[paste0("Row_", i1[ni], "_", as.integer(i2[ni]))]] = pre_prob
    save(ppi_predictProb, file = "../Rdata/14_00_ppi_predictProb.Rdata")
  }
  # indexing
  
  ######################################################
  ###### some PPIs in training and test set are lost because when transform 
  ###### from entrenz id to gene name use different bioMart version
  ###### This step was to trace those PPIs back from training and test set
  
  trainPPI =  get(load('../Rdata/14_00_trainProPair_In1202.Rdata')) 
  valPPI = get(load('../Rdata/14_00_valProPair_In1202.Rdata')) 
  testPPI = get(load('../Rdata/14_00_testProPair_In1202.Rdata'))
  
  #### predicted ppi(losted) in training set 
  index1 = !trainPPI %in% paste0(humapFeature$geneid1, '_', humapFeature$geneid2)
  Training_lostPPI = result$pre_prob_trainSet[index1]
  lostPPI_train = as.data.frame(strSplit(trainPPI[index1], split = '_'), stringsAsFactors = F)
  ## transform entrenz id to gene name
  geneid = unique(c(strSplit(trainPPI[index1], split = '_')[, 1], strSplit(trainPPI[index1], split = '_')[, 2]))
  tmp = ENTREZtoGeneName(geneid)
  tmp = tmp[!duplicated(tmp$entrezgene),]
  rownames(tmp) = tmp$entrezgene
  
  ## 
  tmp1 = tmp[lostPPI_train$V1, ]
  tmp2 = tmp[lostPPI_train$V2, ]
  
  lostPPI_train$geneSymble1 = tmp1$external_gene_name
  lostPPI_train$geneSymble2 = tmp2$external_gene_name
  rownames(lostPPI_train) = trainPPI[index1]
  lostPPI_train$pre_prob = Training_lostPPI
  lostPPI_train = lostPPI_train[, -(1:2)]
  lostPPI_train = lostPPI_train[!is.na(lostPPI_train$geneSymble1) & !is.na(lostPPI_train$geneSymble2), ]
  
  ###### predicted ppi(losted) in test(combined) set
  
  index2 = !c(valPPI, testPPI) %in% paste0(humapFeature$geneid1, '_', humapFeature$geneid2)
  Test_lostPPI = result$pre_prob[index2]
  lostPPI_test = as.data.frame(strSplit(c(valPPI, testPPI)[index2], split = '_'), stringsAsFactors = F)
  ## transform entrenz id to gene name
  geneid = unique(c(lostPPI_test[, 1], lostPPI_test[, 2]))
  tmp = ENTREZtoGeneName(geneid)
  tmp = tmp[!duplicated(tmp$entrezgene),]
  rownames(tmp) = tmp$entrezgene
  
  ## 
  tmp1 = tmp[lostPPI_test$V1, ]
  tmp2 = tmp[lostPPI_test$V2, ]
  
  lostPPI_test$geneSymble1 = tmp1$external_gene_name
  lostPPI_test$geneSymble2 = tmp2$external_gene_name
  rownames(lostPPI_test) = c(valPPI, testPPI)[index2]
  lostPPI_test$pre_prob = Test_lostPPI
  lostPPI_test = lostPPI_test[, -(1:2)]
  lostPPI_test = lostPPI_test[!is.na(lostPPI_test$geneSymble1) & !is.na(lostPPI_test$geneSymble2) & !is.na(lostPPI_test$pre_prob), ]
  
  ###### combine predicted ppi and losted ppi
  load("../Rdata/14_00_ppi_predictProb.Rdata")
  ppi_predictProb = do.call('rbind', ppi_predictProb)
  ppi_predictProb = ppi_predictProb[!is.na(ppi_predictProb$pre_prob), ]
  ppi_predictProb = rbind(ppi_predictProb, lostPPI_train, lostPPI_test)
  save(ppi_predictProb, file = '../Rdata/14_00_ppi_predictProb_all.Rdata')
}


######### split ppis with different top f percentage
# percent = c(0.01, 0.015, 0.02, 0.025, 0.03,0.035, 0.04, 0.05, 0.06)
percent = seq(0.01, 0.2, by = 0.005)
percent = round(dim(ppi_predictProb)[1] * percent)
ppi_predictProb = sortDataframe(ppi_predictProb, by = 'pre_prob', decreasing = T)

for(ni in percent){
  dat = ppi_predictProb[1:ni, ]
  dat0 = data.frame(geneSymble1 = "pro1", geneSymble2 = "pro2", pre_prob = 0)
  dat = rbind(dat0, dat)
  
  write.table(dat, file = paste0('../data/', '14_00_1_', ni, '.tsv'), sep = '\t', 
              row.names = F, col.names = F)
}



############################ Check predict ppis of protein pairs label as 1 (true ppi) ############################ 

## test combine labels
load('../Rdata/12_03_MSBdataExpr_Combine_TextValidation.Rdata')
## train lables
load("../Rdata/10_00_MSBdataExpr_MSBTrain.Rdata") 
## predicted ppi result
result = get(load('../Rdata/ModelSave/12_03/parameter/12_02_model_18_seed186720_350_0.438_140_0.214_25_0.037_Parameter.Rdata'))

## test set predicted ppi 
index = test_combineLables == '1'
ppi_predictTestSet = result$pre_prob[index]
p = hist(ppi_predictTestSet)
plotSave('../figs/14_00_ppi_predictTestSet.png', Plot = p, width = 5, height = 4.5)


## train set predicted ppi 
index = trainlables == '1'
ppi_predictTrainSet = result$pre_prob_trainSet[index]
p = hist(ppi_predictTrainSet)
plotSave('../figs/14_00_ppi_predictTrainingSet.png', Plot = p, width = 5, height = 4.5)

##

















# ## load best model prediction PPI score
# ## best model DL + AP-MS + MS/MS (version 12_03)
# # (12_02_model_18_seed186720_350_0.438_140_0.214_25_0.037_Parameter.Rdata)
# 
# result = get(load('../Rdata/ModelSave/12_03/parameter/12_02_model_18_seed186720_350_0.438_140_0.214_25_0.037_Parameter.Rdata'))
# 
# ## load train PPI 
# train = get(load('../Rdata/14_00_trainProPair_In1202.Rdata')) # in rscript '10_00_CombineMSBdataProteomicData_MSBSets.R'
# trainPPI = data.frame(pro1 = strSplit(train, split = '_')[,1],
#                       pro2 = strSplit(train, split = '_')[,2])
# rownames(trainPPI) = train
# ## match prediction score
# trainPPI$score = result$pre_prob_trainSet
# 
# 
# ## load validation PPI
# val = get(load('../Rdata/14_00_valProPair_In1202.Rdata'))# in rscript '10_00_CombineMSBdataProteomicData_MSBSets.R'
# valPPI = data.frame(pro1 = strSplit(val, split = '_')[,1],
#                     pro2 = strSplit(val, split = '_')[,2])
# rownames(valPPI) = val
# valPPI$score = result$pre_prob
# ##
# PPItotal = rbind(trainPPI, valPPI)
# dim(PPItotal[PPItotal$score > 0.6, ])
# featurePPI =  fread("../rawData_localPC/HuMap_datasets/feature_matrix.txt", 
#                     header = T, select = c(246:247, 245, 248, 249, 265, 277))
# 
# 

