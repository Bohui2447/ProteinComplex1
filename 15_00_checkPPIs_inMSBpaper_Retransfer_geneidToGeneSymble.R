source('00_00_function_complex.R')
options(max.print = 50)
library(funcTools)
library(limma)

{## our assumption is 
  # 1)the MSB paper shows PPI numbers in the paper which generatefrom CORUM complexes,
  # 2)PPIs in feature_matrix.txt are ppis from three AP-MS paper, these ppis are containing AP-MS features 
  # 3) ppis in training and test set are the intersection of ppis in 1 and ppis in 2
  ## this part confirm that our assumption is right
  
  ## total ppis in MSB paper
  TestNeg = read.csv('../rawData_localPC/HuMap_datasets/test_neg_ppis.txt', sep = ' ', header = F)
  TestPos = read.csv('../rawData_localPC/HuMap_datasets/test_positive_ppis.txt', sep = '\t', header = F)
  TrainNeg = read.csv('../rawData_localPC/HuMap_datasets/train_neg_ppis.txt', sep = ' ', header = F)
  trainPos = read.csv('../rawData_localPC/HuMap_datasets/train_positive_ppis.txt', sep = ' ', header = F)
  
  totalPPI = rbind(TestNeg, TestPos, TrainNeg,trainPos)
  totalPPI$pro1pro2 = paste0(totalPPI$V1, '_', totalPPI$V2)
  totalPPI$pro2pro1 = paste0(totalPPI$V2, '_', totalPPI$V1)
  
  totalPPIpro1pro2 = unique(c(totalPPI$pro1pro2, totalPPI$pro2pro1))
  
  ## PPIs with AP-MS feature in MSB paper
  featurePPI =  fread("../rawData_localPC/HuMap_datasets/feature_matrix.txt", 
                      header = T, select = c(246:247))
  featurePPI =  fread("../rawData_localPC/HuMap_datasets/feature_matrix.txt", 
                      header = T, select = c(246:247, 245, 248, 249, 265, 277))
  
  tmp = featurePPI$frozenset_geneids_str_order
  tmp2 = strSplit(gsub("'|,|\\[|\\]", "", tmp), " ")
  tmp3 = t(apply(tmp2, 1, order))
  
  tmp4 = t(apply(featurePPI[,1:2], 1, order))
  table(tmp4[,1])
  table(tmp4[,2])
  
  featurePPI = as.data.frame(featurePPI)
  featurePPI[,1] = as.integer(featurePPI[,1])
  featurePPI[,2] = as.integer(featurePPI[,2])
  head(featurePPI, 20)
  
  featurePPI$pro1pro2 = paste0(featurePPI[,1], '_', featurePPI[,2])
  tmp = intersect(featurePPI$pro1pro2, totalPPIpro1pro2)
  length(tmp)
}


######### generate new protein pair for PPI prediction
if(file.exists('../Rdata/15_00_humapFeature_transformGeneidToGenesymble.Rdata')){
  cat('load ../Rdata/15_00_humapFeature_transformGeneidToGenesymble.Rdata')
  load('../Rdata/15_00_humapFeature_transformGeneidToGenesymble.Rdata')
}else{
  # AP-MS data (MSB paper data)
  HumapFeature0 = get(load("../RdataAll/11featureCleanHuMap.Rdata"))
  
  # grep ENSG id in entrenz id (geneid1 and geneid2)
  ENSG1 = HumapFeature0[grep("^[A-z]+", HumapFeature0$geneid1, value = F),]$geneid1
  ENSG2 = HumapFeature0[grep("^[A-z]+", HumapFeature0$geneid2, value = F),]$geneid2
  
  # transform ENSG to entrenzgene id
  ENSG1Entrenzgene_89 = ENSGtoENTREZ(ENSG1)
  ENSG1Entrenzgene_89 = ENSG1Entrenzgene_89[!is.na(ENSG1Entrenzgene_89$entrezgene),]
  ENSG1Entrenzgene_89 = ENSG1Entrenzgene_89[!duplicated(ENSG1Entrenzgene_89$ensembl_gene_id),]
  
  
  ENSG2Entrenzgene_89 = ENSGtoENTREZ(unique(ENSG2))
  ENSG2Entrenzgene_89 = ENSG2Entrenzgene_89[!is.na(ENSG2Entrenzgene_89$entrezgene), ]
  ENSG2Entrenzgene_89 = ENSG2Entrenzgene_89[!duplicated(ENSG2Entrenzgene_89$ensembl_gene_id),]
  
  ENSGEntrenzgene_89 = rbind(ENSG1Entrenzgene_89, ENSG2Entrenzgene_89)
  ENSGEntrenzgene_89 = ENSGEntrenzgene_89[!duplicated(ENSGEntrenzgene_89),]
  
  ## split humapfeature0
  ensgRows = sort(union(grep("^[A-z]+", HumapFeature0$geneid1), 
                        grep("^[A-z]+", HumapFeature0$geneid2)))
  HumapFeature1 = HumapFeature0[ensgRows,]
  
  # match  
  rownames(ENSGEntrenzgene_89) = ENSGEntrenzgene_89$ensembl_gene_id
  tmp1 = ENSGEntrenzgene_89[HumapFeature1$geneid1,]
  tmp2 = ENSGEntrenzgene_89[HumapFeature1$geneid2,]
  
  # replace ensembl gene id with entrenzgene  in humapfeature1
  ensgindex1 = grep('ENSG', HumapFeature1$geneid1)
  HumapFeature1[ensgindex1, ]$geneid1 = tmp1[ensgindex1, ]$entrezgene
  
  ensgindex2 = grep('ENSG', HumapFeature1$geneid2)
  HumapFeature1[ensgindex2, ]$geneid2 = tmp2[ensgindex2, ]$entrezgene
  ## remove rows in geneid with NA
  HumapFeature1 = HumapFeature1[!(is.na(HumapFeature1$geneid1) | is.na(HumapFeature1$geneid2)), ]
  ## sort geneid in humapfeature1
  tmp = HumapFeature1[, 1:2]
  tmp1 = t(apply(tmp, 1, sort))
  HumapFeature1$geneid1 = unname(tmp1[,1])
  HumapFeature1$geneid2 = unname(tmp1[,2])
  
  ## combine humapfeature0 (remove rows with ENSG id) and humapfeagure1
  humapFeature = rbind(HumapFeature0[-ensgRows, ], HumapFeature1)
  
  ##### transform entrenzgene id to gene symble in humapFeature
  geneid = unique(c(unique(humapFeature$geneid1), unique(humapFeature$geneid2)))
  
  geneid_genename = ENTREZtoGeneName(geneid)
  # remove duplicated geneid and gene symble
  geneid_genename = geneid_genename[!duplicated(geneid_genename$entrezgene),]
  geneid_genename = geneid_genename[!duplicated(geneid_genename$external_gene_name),]
  
  rownames(geneid_genename) = geneid_genename$entrezgene
  
  tmp1 = geneid_genename[humapFeature$geneid1, ]
  tmp2 = geneid_genename[humapFeature$geneid2, ]
  
  humapFeature$geneSymble1 = tmp1$external_gene_name
  humapFeature$geneSymble2 = tmp2$external_gene_name
  
  # remove humapFeature gene symble contain NA
  humapFeature = humapFeature[!(is.na(humapFeature$geneSymble1) | is.na(humapFeature$geneSymble2)), ]
  # remove duplicated ppis
  xx = paste0(humapFeature$geneSymble1, '_', humapFeature$geneSymble2)
  humapFeature$pro1_pro2 = xx
  humapFeature = humapFeature[!duplicated(humapFeature$pro1_pro2),]
  humapFeature = humapFeature[, -267]
  
  save(humapFeature, file = '../Rdata/15_00_humapFeature_transformGeneidToGenesymble.Rdata')
}





# 
# featurePPI =  fread("../rawData_localPC/HuMap_datasets/feature_matrix.txt", 
#                     header = T, nrows = 100)








