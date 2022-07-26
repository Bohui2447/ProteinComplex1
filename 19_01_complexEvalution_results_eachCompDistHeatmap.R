library(funcTools)
library(ggplot2)
library(pheatmap)
options(max.print = 50)

## plot and save complex gene order in heatmap 
# load PXD_combine dataset and complexes
PXD_combine = get(load("../Rdata_localPC/3PXD_combine.Rdata"))
# log data
logExpr = log10(PXD_combine + 1)
# load predicted complex
# fl = c('14_00_1_167384.tsv_0.2_0.8_4_MCL.txt', '14_00_1_131516.tsv_0.25_0.8_5_MCL.txt',  '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt')
comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', '14_00_1_131516.tsv_0.25_0.8_5_MCL.txt'), 
                header = F, sep = '\t', stringsAsFactors = F)
# load distance matrix 
load('../Rdata/19_00_14_00_1_131516.tsv_0.25_0.8_5_MCL.txt_distance.Rdata')
distance = as.matrix(di)
proOrder = list()
for(ni in 1:nrow(comp)){
  comp_ni = comp[ni,]
  pro_ni = t(strSplit(comp_ni, split = ' '))[,1]
  # distance plot
  compDist = distance[pro_ni, pro_ni]
  
  p = pheatmap(compDist, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
               border_color = NA,  cluster_cols = T, cluster_rows = T)
  plotSave(filename = paste0('../figs/19_01_comp_14_00_1_131516.tsv_0.25_0.8_5_MCL.txt/complex_', ni, '.png'), 
           Plot = p,width = 6, height = 6, dpi = 300)
  proOrder[[ni]] =  rownames(compDist)[p$tree_row$order]
  
}
save(proOrder, file = '../Rdata/19_01_eachCompGeneOrder_14_00_1_131516.tsv_0.25_0.8_5_MCL.txt.Rdata')


c('COX7A2', 'COX7A2L', 'MAGEA2', 'MAGEA2B', 'TOMM70', 'TOMM70A', 'MYEOV', 'MYEOV2', 'MRE11', 'MRE11A')
pro_ni[!pro_ni %in% rownames(distance)]
pro_ni[!pro_ni %in% rownames(distance)] = 'MRE11A'



# fl = c('14_00_1_167384.tsv_0.2_0.8_4_MCL.txt', '14_00_1_131516.tsv_0.25_0.8_5_MCL.txt',  '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt')
comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', '14_00_1_167384.tsv_0.2_0.8_4_MCL.txt'), 
                header = F, sep = '\t', stringsAsFactors = F)
proAll = paste(comp[,1], collapse = ' ')
proAll = t(strSplit(proAll, split = ' '))
proAll = unique(proAll[,1])
# load distance matrix 
load('../Rdata/19_00_14_00_1_167384.tsv_0.2_0.8_4_MCL.txt_distance.Rdata')
distance = as.matrix(di)
## replace gene name with alternative name in distance matrix
if(FALSE){
# proAll[!(proAll %in% rownames(distance) )] 
# c("GCN1" "RPS17" "MAGEA2" "UTP11"  "UFD1"   "CRAMP1" "COX7A2" "FAM92A" "TOMM70" "MYEOV"  "MRE11")s
rownames(distance)[!(rownames(distance) %in% proAll)] = c("GCN1", "RPS17", "MAGEA2" ,"UTP11" ,
                                                          "UFD1"  , "CRAMP1" ,"COX7A2", "FAM92A", "TOMM70" ,"MYEOV" , "MRE11")
colnames(distance)[!(colnames(distance) %in% proAll)] = c("GCN1", "RPS17", "MAGEA2" ,"UTP11" ,
                                                          "UFD1"  , "CRAMP1" ,"COX7A2", "FAM92A", "TOMM70" ,"MYEOV" , "MRE11")
# c("GCN1L1"    "RPS17L"    "MAGEA2B"   "UTP11L"    "UFD1L"     "CRAMP1L" "COX7A2L.1"  "FAM92A1"    "TOMM70A"   "MYEOV2"    "MRE11A"  )
}
proOrder = list()
for(ni in 1:nrow(comp)){
  comp_ni = comp[ni,]
  pro_ni = t(strSplit(comp_ni, split = ' '))[,1]
  # distance plot
  compDist = distance[pro_ni, pro_ni]
  
  p = pheatmap(compDist, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
               border_color = NA,  cluster_cols = T, cluster_rows = T)
  plotSave(filename = paste0('../figs/19_01_comp_14_00_1_167384.tsv_0.2_0.8_4_MCL.txt/complex_', ni, '.png'), 
           Plot = p,width = 6, height = 6, dpi = 300)
  proOrder[[ni]] =  rownames(compDist)[p$tree_row$order]
  
}
save(proOrder, file = '../Rdata/19_01_eachCompGeneOrder_14_00_1_167384.tsv_0.2_0.8_4_MCL.txt.Rdata')




## 14_00_1_119560.tsv_0.2_0.8_5_MCL.txt
comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt'), 
                header = F, sep = '\t', stringsAsFactors = F)
proAll = paste(comp[,1], collapse = ' ')
proAll = t(strSplit(proAll, split = ' '))
proAll = unique(proAll[,1])
# load distance matrix 
load('../Rdata/19_00_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt_distance.Rdata')
distance = as.matrix(di)
## replace gene name with alternative name in distance matrix
if(FALSE){
  # proAll[!(proAll %in% rownames(distance) )] 
  # c( "GCN1"   "RPS17"  "MAGEA2" "UTP11"  "UFD1"   "TOMM70" "CRAMP1" "COX7A2" "MRE11"  "FAM92A" "MYEOV" )
  rownames(distance)[!(rownames(distance) %in% proAll)] = c( "GCN1" , "RPS17","MAGEA2","UTP11","UFD1" , "TOMM70","CRAMP1","MRE11","FAM92A","COX7A2","MYEOV")
  colnames(distance)[!(colnames(distance) %in% proAll)] = c( "GCN1" , "RPS17","MAGEA2","UTP11","UFD1" , "TOMM70","CRAMP1","MRE11","FAM92A","COX7A2","MYEOV")
  # c("GCN1L1"    "RPS17L"    "MAGEA2B"   "UTP11L"    "UFD1L"     "TOMM70A"   "CRAMP1L"   "MRE11A"    "FAM92A1"   "COX7A2L.1" "MYEOV2"   )
}
proOrder = list()
for(ni in 1:nrow(comp)){
  comp_ni = comp[ni,]
  pro_ni = t(strSplit(comp_ni, split = ' '))[,1]
  # distance plot
  compDist = distance[pro_ni, pro_ni]
  
  p = pheatmap(compDist, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
               border_color = NA,  cluster_cols = T, cluster_rows = T)
  plotSave(filename = paste0('../figs/19_01_comp_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt/complex_', ni, '.png'), 
           Plot = p,width = 6, height = 6, dpi = 300)
  proOrder[[ni]] =  rownames(compDist)[p$tree_row$order]
  
}
save(proOrder, file = '../Rdata/19_01_eachCompGeneOrder_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt.Rdata')




########## protein numbers in 3 predict complex dataset

fl = c('14_00_1_167384.tsv_0.2_0.8_4_MCL.txt', '14_00_1_131516.tsv_0.25_0.8_5_MCL.txt',  '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt')
for(mi in fl){
  comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', mi), 
                  header = F, sep = '\t', stringsAsFactors = F)
  proNumPerComp = apply(comp, 1, function(x) length(t(strSplit(x, split = ' '))[,1]))
  {
    xxx = as.data.frame(table(proNumPerComp))
    xxx$proNumPerComp = as.numeric(as.character(xxx$proNumPerComp))
    `25-30` = data.frame(proNumPerComp = "25-30", 
                         Freq = sum(subset(xxx, proNumPerComp >= 25 & proNumPerComp <= 30)$Freq))
    `31-35` = data.frame(proNumPerComp = "31-35", 
                         Freq = sum(subset(xxx, proNumPerComp >= 31 & proNumPerComp <= 35)$Freq))
    `>35` = data.frame(proNumPerComp = ">35", 
                       Freq = sum(subset(xxx, proNumPerComp > 35)$Freq))
    freqComplex = rbind(subset(xxx, proNumPerComp < 25),
                        `25-30`, `31-35`, `>35`)
    freqComplex$proNumPerComp = factor(freqComplex$proNumPerComp,
                                       freqComplex$proNumPerComp)
  }
    bks = c(seq(1,10,by=1), seq(10,100,by=10), 
            seq(100,1000,by=100), seq(1000,10000,by=1000))
    lbs = bks
    lbs[!bks %in% c(1, 10, 100, 1000, 2000, 3000)] = ""
    p = ggplot(freqComplex, aes(proNumPerComp, Freq)) + 
      geom_bar(stat = "identity") + 
      theme_Publication(x_angle = 90, x_hjust = 1) +
      scale_y_log10(breaks= bks, labels = lbs)+
      ylab('Frequency') + xlab('# Proteins') +
      geom_text(aes(label = Freq, angle = 90, hjust = 1.3), color = "yellow")
    plotSave(filename = paste0('../figs/19_01_proNumPerComp_', mi, '.png'), 
             Plot = p,width = 6, height = 6, dpi = 300)
    plotSave(filename = paste0('../figs/19_01_proNumPerComp_', mi, '.pdf'), 
             Plot = p,width = 6, height = 6, dpi = 300)
  
}



########## transpose complexes in column
fl = c('14_00_1_167384.tsv_0.2_0.8_4_MCL.txt', '14_00_1_131516.tsv_0.25_0.8_5_MCL.txt',  '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt')
for(mi in fl){
  comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', mi), 
                  header = F, sep = '\t', stringsAsFactors = F)
  proNumPerComp = apply(comp, 1, function(x) length(t(strSplit(x, split = ' '))[,1]))
  tmp = matrix(data = NA, nrow = max(proNumPerComp), ncol = nrow(comp))
  for(ni in 1:nrow(comp)){
    comp_ni = comp[ni,]
    pro_ni = t(strSplit(comp_ni, split = ' '))[,1]
    tmp[,ni] = c(pro_ni, rep(NA, nrow(tmp) - length(pro_ni)))
  }
  write.csv(tmp, file = paste0('../Rdata/19_00_complexEvalution_result/complexInColumn_', mi, '.csv'))
}




########## plot select complexes in 14_00_1_119560.tsv_0.2_0.8_5_MCL.txt
comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt'), 
                header = F, sep = '\t', stringsAsFactors = F)
proAll = paste(comp[,1], collapse = ' ')
proAll = t(strSplit(proAll, split = ' '))
proAll = unique(proAll[,1])
# load distance matrix 
load('../Rdata/19_00_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt_distance.Rdata')
distance = as.matrix(di)
## replace gene name with alternative name in distance matrix
if(FALSE){
  # proAll[!(proAll %in% rownames(distance) )] 
  # c( "GCN1"   "RPS17"  "MAGEA2" "UTP11"  "UFD1"   "TOMM70" "CRAMP1" "COX7A2" "MRE11"  "FAM92A" "MYEOV" )
  rownames(distance)[!(rownames(distance) %in% proAll)] = c( "GCN1" , "RPS17","MAGEA2","UTP11","UFD1" , "TOMM70","CRAMP1","MRE11","FAM92A","COX7A2","MYEOV")
  colnames(distance)[!(colnames(distance) %in% proAll)] = c( "GCN1" , "RPS17","MAGEA2","UTP11","UFD1" , "TOMM70","CRAMP1","MRE11","FAM92A","COX7A2","MYEOV")
  # c("GCN1L1"    "RPS17L"    "MAGEA2B"   "UTP11L"    "UFD1L"     "TOMM70A"   "CRAMP1L"   "MRE11A"    "FAM92A1"   "COX7A2L.1" "MYEOV2"   )
}
# load protein order in every complex
proOrder = get(load('../Rdata/19_01_eachCompGeneOrder_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt.Rdata'))

# comp57 = t(strSplit(comp[57, ], split = ' '))[,1]
# pro = c(proOrder[[3195]], proOrder[[571]], proOrder[[1404]], proOrder[[76]])
comp_potential = c(76,  15, 124, 163, 131, 401, 509, 502, 571, 748, 739, 1162, 1231, 1404, 1688, 3195)
compCombn = t(combn(comp_potential, m = 4))
for(mi in 1:nrow(compCombn)){
comb_mi = compCombn[mi, ]
 
  pro = c(proOrder[[comb_mi[1]]], proOrder[[comb_mi[2]]], proOrder[[comb_mi[3]]], proOrder[[comb_mi[4]]])
# distance plot
compDist = distance[pro, pro]

# generate dataframe for annotation
annoRow = data.frame(Complex = c(rep('Complex1', length(proOrder[[comb_mi[1]]])),
                                 rep('Complex2', length(proOrder[[comb_mi[2]]])),
                                 rep('Complex3', length(proOrder[[comb_mi[3]]])),
                                 rep('complex4', length(proOrder[[comb_mi[4]]]))))
rownames(annoRow) = pro
# p = pheatmap(compDist, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
             # border_color = NA, annotation_row = annoRow, cluster_cols = T, cluster_rows = T)
p = pheatmap(compDist, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
             border_color = NA, annotation_row = annoRow, cluster_cols = F, cluster_rows = F)
out = paste(comb_mi[1], comb_mi[2], comb_mi[3], comb_mi[4], sep = '_')
plotSave(paste0("../figs/19_01_compDistanceHeatmap/comp_", out, '.png'), Plot = p, width = 7, height = 7)
}





########## plot select complexes in 14_00_1_119560.tsv_0.2_0.8_5_MCL.txt (#### Used in paper)
comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt'), 
                header = F, sep = '\t', stringsAsFactors = F)
proAll = paste(comp[,1], collapse = ' ')
proAll = t(strSplit(proAll, split = ' '))
proAll = unique(proAll[,1])
# load distance matrix 
load('../Rdata/19_00_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt_distance.Rdata')
distance = as.matrix(di)
## replace gene name with alternative name in distance matrix
if(FALSE){
  # proAll[!(proAll %in% rownames(distance) )] 
  # c( "GCN1"   "RPS17"  "MAGEA2" "UTP11"  "UFD1"   "TOMM70" "CRAMP1" "COX7A2" "MRE11"  "FAM92A" "MYEOV" )
  rownames(distance)[!(rownames(distance) %in% proAll)] = c( "GCN1" , "RPS17","MAGEA2","UTP11","UFD1" , "TOMM70","CRAMP1","MRE11","FAM92A","COX7A2","MYEOV")
  colnames(distance)[!(colnames(distance) %in% proAll)] = c( "GCN1" , "RPS17","MAGEA2","UTP11","UFD1" , "TOMM70","CRAMP1","MRE11","FAM92A","COX7A2","MYEOV")
  # c("GCN1L1"    "RPS17L"    "MAGEA2B"   "UTP11L"    "UFD1L"     "TOMM70A"   "CRAMP1L"   "MRE11A"    "FAM92A1"   "COX7A2L.1" "MYEOV2"   )
}
# load protein order in every complex
proOrder = get(load('../Rdata/19_01_eachCompGeneOrder_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt.Rdata'))
comp_potential = c(76, 57, 15, 124, 163, 131, 401, 509, 502, 571, 748, 739, 1162, 1231, 
                   1404, 1406, 1569, 1688, 1862, 1866,3195, 4407, 4510, 4689)

for(mi in comp_potential){
  pro = (proOrder[[mi]])
  # distance plot
  compDist = distance[pro, pro]
  p = pheatmap(compDist, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
               border_color = NA,cluster_cols = T, cluster_rows = T, 
               show_rownames = T, show_colnames = T)
  
  plotSave(paste0("../figs/19_01_1_119560_comp_", mi, '.png'), Plot = p, width = 7, height = 7)
  plotSave(paste0("../figs/19_01_1_119560_comp_", mi, '.pdf'), Plot = p, width = 7, height = 7)
}













