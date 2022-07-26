# 'Rdata_localPC' is 'RdataAll'

options(max.print = 100)

if(file.exists("../Rdata_localPC/5exprLabel_ratio.Rdata")){
  cat("load ../Rdata_localPC/5exprLabel_ratio.Rdata")
  load("../Rdata_localPC/5exprLabel_ratio.Rdata")
} else {
  #### load PXD_combine dataset and complexes
  PXD_combine = get(load("../Rdata_localPC/3PXD_combine.Rdata")) ## new dataset
  comProPro = get(load("../Rdata_localPC/4comProPro_CORUM.Rdata"))
  
  
  # log data
  logExpr = log10(PXD_combine + 1)
  # tmp_na = apply(logExpr, 1, function(x) sum(is.na(x)))
  
  logExpr1 = logExpr[comProPro$protein1,]
  logExpr2 = logExpr[comProPro$protein2,]
  
  # test na in dataset
  tmp_na = apply(logExpr1, 1, function(x) sum(is.na(x))) ## proteins not in complexes will be NA
  
  rownames(logExpr1) = paste0(comProPro$protein1, "_", comProPro$protein2)
  rownames(logExpr2) = paste0(comProPro$protein1, "_", comProPro$protein2)
  
  rowNa = (is.na(logExpr1[,1]) | is.na(logExpr2[,1]))
  
  # exprAbs = abs(logExpr1[!rowNa,] - logExpr2[!rowNa,]) ## change to eucledian distance, to be done ***
  exprAbs = abs((logExpr1[!rowNa,] - logExpr2[!rowNa,])/(logExpr1[!rowNa,] + logExpr2[!rowNa,] + 1))
  
  
  # total proteins from complex
  comProTotal = unique(c(comProPro$protein1, comProPro$protein2))
  
  #### generate negative protein pairs
  set.seed(1234) # no seed was set in the initial codes, so the data might change
  randomNeg = t(sapply(1:(3*nrow(logExpr1)), function(x){
    sample(comProTotal, size = 2, replace = FALSE)}
  ))
  
  # remove duplicated combinations
  neg = t(apply(randomNeg, 1, sort))
  neg = neg[!duplicated(neg),]
  
  # remove the combinations in positive set
  neg = neg[!(paste0(neg[,1], "_", neg[,2]) %in% rownames(logExpr1)), ]
  
  # the number of negative set is the same as that of positive set
  negSet = neg[1:(nrow(logExpr1)), ]
  
  ###### negative dataset ######
  
  logExpr1_neg = logExpr[negSet[,1],]
  logExpr2_neg = logExpr[negSet[,2],]
  
  rownames(logExpr1_neg) = rownames(logExpr2_neg) = paste0(negSet[,1], "_", negSet[,2])
  
  rowNa = (is.na(logExpr1_neg[,1]) | is.na(logExpr2_neg[,1]))
  
  # exprAbs_neg = abs(logExpr1_neg[!rowNa,] - logExpr2_neg[!rowNa,]) # change to abs((A-B)/(A+B+1))
  exprAbs_neg =   abs((logExpr1_neg[!rowNa,] - logExpr2_neg[!rowNa,])/
                        (logExpr1_neg[!rowNa,] + logExpr2_neg[!rowNa,] + 1))
  
  
  ##### combine positive dataset and negative dataset #####
  expr = rbind(exprAbs, exprAbs_neg)
  colnames(expr) = paste0("V", 1:ncol(expr))
  exprLabel = cbind(expr, label = c(rep(1, nrow(exprAbs)), rep(0, nrow(exprAbs_neg))))
  exprLabel = as.data.frame(exprLabel)
  
  ## delete all 0 columns
  # index = read.csv("../Rdata_localPC/deleteCol.csv", header = T, stringsAsFactors = F, as.is = T)
  # exprLabel = exprLabel[, -index$delete1]
  # write.csv(exprLabel, file = "../1rawData/exprLabel.csv")
  save(exprLabel, file = "../Rdata_localPC/5exprLabel_ratio.Rdata") ### this will be different from the original data, because the random seed was not set.
}