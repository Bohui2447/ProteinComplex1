library(funcTools)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(pracma)
options(max.print = 50)

###### combine all complex evalution result by k-clique
if(file.exists('../Rdata/19_00_complexEvalution_result/compEva_result.csv')){
  cat('read.csv(../Rdata/19_00_complexEvalution_result/compEva_result.csv)')
  read.csv('../Rdata/19_00_complexEvalution_result/compEva_result.csv', header = T, row.names = 1)
} else{
  path = "../Rdata/17_00_MCL_onClusterOneResult_txt/"
  file = dir(path, "*output.tsv")
  
  compEva_result = list()
  for(mi in file){
    fileName = paste0(path, mi)
    res = read.csv(file = fileName, header = T, sep = "\t")
    compEva_result[[mi]] = res
    
  }
  compEva_result = do.call(rbind, compEva_result)
  ## save complex evaluation result
  Folder = '../Rdata/19_00_complexEvalution_result'
  if(!dir.exists(Folder)){
    dir.create(Folder)
  }
  write.csv(compEva_result, file = '../Rdata/19_00_complexEvalution_result/compEva_result.csv')
}


########################  plot using F-grand k-clique
## read.csv(file = fileName, header = T, sep = "\t") '14_00_1_442373.tsv_0.2_0.8_5_MCL.txt_output.tsv'
compEva_result = read.csv('../Rdata/19_00_complexEvalution_result/compEva_result.csv', header = T, row.names = 1)
tmp = strSplit(rownames(compEva_result), split = '_')

NumPPI = as.numeric(strSplit(tmp[,4], split = '\\.')[,1])
minDensity = as.numeric(tmp[,5])
maxOverlap = tmp[,6]
inflation = as.numeric(tmp[,7])
datPlot = data.frame(NumPPI = NumPPI, minDensity = minDensity,
                     maxOverlap = maxOverlap, inflation = inflation,
                     Fgrand = compEva_result$F.Grand.K.Clique, 
                     param = paste(NumPPI, minDensity, maxOverlap, inflation, sep = '_'), 
                     stringsAsFactors = F)
datPlot = sortDataframe(datPlot, by = 'Fgrand', decreasing = F)
datPlot$Index = order(datPlot$Fgrand, decreasing = F)
# tmp = datPlot[datPlot$NumPPI == 442373, ]

# ggplot(datPlot, aes(x = NumPPI, y = Fgrand)) +
#   geom_point(aes( shape = maxOverlap, size = minDensity , color = inflation), alpha = 0.2) # 
p = ggplot(datPlot, aes(x = Index, y = Fgrand)) +
  geom_point(data = datPlot[datPlot$Index == 7598,], 
             aes(x = Index, y = Fgrand),shape=23, fill="red", color="darkred", size=4)+
  geom_point(alpha = 0.05)+
  theme_Publication()
print(p)

ggsave(p, filename = '../figs/19_00_F_grand.png', width = 6, height = 5, dpi = 300)
ggsave(p, filename = '../figs/19_00_F_grand.pdf', width = 6, height = 5, dpi = 300)



################## Check complexes number 
if(file.exists('../Rdata/19_00_complexEvalution_result/compEva_result_complexNumber.csv')){
  cat('read.csv(../Rdata/19_00_complexEvalution_result/compEva_result_complexNumber.csv)')
  read.csv('../Rdata/19_00_complexEvalution_result/compEva_result_complexNumber.csv', header = T, row.names = 1)
}else {
  fl = dir('../Rdata/17_00_MCL_onClusterOneResult_txt/', pattern = 'MCL.txt')
  fl = fl[seq(1,length(fl), by = 2)]
  compNum = list()
  for(mi in fl){
    comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', mi), 
                    header = F, sep = '\t')
    numberComp = nrow(comp)
    proAll = paste(comp[,1], collapse = ' ')
    proAll = t(strSplit(proAll, split = ' '))
    numberPro = length(unique(proAll[,1]))
    compNum[[mi]] = c(numberComp = numberComp,numberPro = numberPro)
    
  }
  compNum = do.call(rbind, compNum)
  rownames(compNum) = paste0(rownames(compNum), '_output.tsv')
  
  compEva_result = read.csv('../Rdata/19_00_complexEvalution_result/compEva_result.csv', header = T, row.names = 1)
  
  tmp = compEva_result[rownames(compNum),]
  tmp$complexNum = compNum[,1]
  tmp$proteinNum = compNum[,2]
  write.csv(tmp, file = '../Rdata/19_00_complexEvalution_result/compEva_result_complexNumber.csv')
}



################## Complexes number F-grand plot
## 
tmp = read.csv('../Rdata/19_00_complexEvalution_result/compEva_result_complexNumber.csv', header = T, row.names = 1)
tmp1 = tmp[, c('F.Grand.K.Clique', 'complexNum','proteinNum', 'Clique.Weighted.hmean..F.weighted.K.Clique.')]
p = ggplot(tmp1, aes(F.Grand.K.Clique, complexNum))+
  geom_point(alpha = 0.1) + #, aes(size = log10(tmp1$proteinNum))
  xlab('F-grand') + ylab('# Complexes')+
  theme_Publication()
print(p)
ggsave(p, filename = '../figs/19_00_Fgrand_ComplexNumber.png', width = 5.5, height = 5, dpi = 300)
ggsave(p, filename = '../figs/19_00_Fgrand_ComplexNumber.pdf', width = 5.5, height = 5, dpi = 300)



################## Complexes number F-grand plot
## 
tmp = read.csv('../Rdata/19_00_complexEvalution_result/compEva_result_complexNumber.csv', header = T, row.names = 1)
tmp1 = tmp[, c('F.Grand.K.Clique', 'complexNum','proteinNum', 'Clique.Weighted.hmean..F.weighted.K.Clique.')]
tmp1$interval = cut(tmp1$complexNum, breaks = seq(0, 6000, by = 1000))
p = ggplot(tmp1, aes(F.Grand.K.Clique, complexNum))+
  geom_point(alpha = 0.1, aes(color = tmp1$interval)) + #, aes(size = log10(tmp1$proteinNum))
  xlab('F-grand') + ylab('# Complexes')+
  theme_Publication()
print(p)
ggsave(p, filename = '../figs/19_00_Fgrand_ComplexNumber_color.png', width = 5.5, height = 5, dpi = 300)
ggsave(p, filename = '../figs/19_00_Fgrand_ComplexNumber_color.pdf', width = 5.5, height = 5, dpi = 300)




################## Plot predicted complexes correlation and distance using MS/MS data(feature)
##  best complex result candidate
# index = rownames(tmp1[tmp1$F.Grand.K.Clique> 0.46 &tmp1$complexNum > 5000, ])
# index = rownames(tmp1[tmp1$F.Grand.K.Clique> 0.46 &tmp1$complexNum > 4600, ])

## CORUM complexes
corum = read.csv(file = '../Rdata/17_00_MCL_onClusterOneResult_txt/10_CorumComp.txt', header = F)


#### caculate and save distance matrix 
# load PXD_combine dataset and complexes
PXD_combine = get(load("../Rdata_localPC/3PXD_combine.Rdata"))
# log data
logExpr = log10(PXD_combine + 1)
if(file.exists('../Rdata/19_00_distance_all.Rdata')){
  load('../Rdata/19_00_distance_all.Rdata')
  cat('load ../Rdata/19_00_distance_all.Rdata')
} else {
  distance = distmat(logExpr, logExpr)
  save(distance, file = '../Rdata/19_00_distance_all.Rdata')
}

fl = c('14_00_1_167384.tsv_0.2_0.8_4_MCL.txt',
       '14_00_1_131516.tsv_0.25_0.8_5_MCL.txt', 
       '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt')
for(mi in fl){
  ## predict complex
  comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', mi), 
                  header = F, sep = '\t', stringsAsFactors = F)
  # caculate distance matrix
  proAll = paste(comp[,1], collapse = ' ')
  proAll = t(strSplit(proAll, split = ' '))
  compExp = logExpr[unique(proAll[,1]), ]
  di = dist(compExp)
  # save distance
  save(di, file =  paste0('../Rdata/19_00_', mi, '_distance.Rdata'))
  # plot distance heatmap 
  distance = as.matrix(di)
  p = pheatmap(distance, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
               border_color = NA, show_rownames = F, show_colnames = F,
               cluster_cols = T, cluster_rows = T)
  save(p, file = paste0('../Rdata/19_00_', mi, '_heatmap.Rdata'))
  plotSave(paste0("../figs/19_00_", mi, "_CompDist_pheatmap.png"), Plot = p, width = 12, height = 12)
  # heatmap genes in order
  dat = data.frame(index = 1:nrow(distance), pro = rownames(distance)[p$tree_row$order])
  save(dat, file = paste0('../Rdata/19_00_', mi, '_heatmapGeneOrder.Rdata'))
}



########### save mean, minimum, maxmum distance for each complex
# predicted complexes
comp = read.csv(file = paste0('../Rdata/17_00_MCL_onClusterOneResult_txt/', '14_00_1_119560.tsv_0.2_0.8_5_MCL.txt'), 
                header = F, sep = '\t', stringsAsFactors = F)
# all protein distance using ms/ms feature
distance = get(load('../Rdata/19_00_14_00_1_119560.tsv_0.2_0.8_5_MCL.txt_distance.Rdata'))
di = as.matrix(distance)

if(FALSE){
  proAll = paste(comp[,1], collapse = ' ')
  proAll = t(strSplit(proAll, split = ' '))
  proAll = unique(proAll[,1])
  # proAll[!(proAll %in% rownames(distance) )] 
  # c("GCN1" "RPS17" "MAGEA2" "UTP11"  "UFD1"   "CRAMP1" "COX7A2" "FAM92A" "TOMM70" "MYEOV"  "MRE11")s
  rownames(di)[!(rownames(di) %in% proAll)] = c("GCN1", "RPS17", "MAGEA2" ,"UTP11" ,
                                                            "UFD1"  , "CRAMP1" ,"COX7A2", "FAM92A", "TOMM70" ,"MYEOV" , "MRE11")
  colnames(di)[!(colnames(di) %in% proAll)] = c("GCN1", "RPS17", "MAGEA2" ,"UTP11" ,
                                                            "UFD1"  , "CRAMP1" ,"COX7A2", "FAM92A", "TOMM70" ,"MYEOV" , "MRE11")
  # c("GCN1L1"    "RPS17L"    "MAGEA2B"   "UTP11L"    "UFD1L"     "CRAMP1L" "COX7A2L.1"  "FAM92A1"    "TOMM70A"   "MYEOV2"    "MRE11A"  )
}

## complex distance
if(file.exists('../Rdata/19_00_complex_distance.Rdata')){
  load('../Rdata/19_00_complex_distance.Rdata')
  cat('load ../Rdata/19_00_complex_distance.Rdata')
}else{
  complex_distance = list()
  for(mi in 1:nrow(comp)){
    pro = strSplit(comp[mi, ], split = ' ')[1,]
    compDist = di[pro, pro]
    upperDist = compDist[upper.tri(compDist)]
    dist_average = mean(upperDist)
    dist_min = min(upperDist)
    dist_max = max(upperDist)
    dist_median = median(upperDist)
    proNum = length(pro)
    complex_distance[[mi]] = c(proNum = proNum, dist_average = dist_average, dist_min = dist_min, 
                               dist_max = dist_max, dist_median = dist_median)
  }
  
  complex_distance = do.call('rbind', complex_distance)
  write.csv(complex_distance, file = '../Rdata/19_00_complex_distance.csv')
  save(complex_distance, file = '../Rdata/19_00_complex_distance.Rdata')
}


## plot
{
dat = data.frame(complex_distance)
datMelt = melt(dat, measure.vars = c('dist_average',   'dist_min',   'dist_max', 'dist_median'))
g = ggplot(data = datMelt, aes(variable, value)) +
  geom_boxplot() + 
  labs(x = '', y = 'Distance', title = 'Complexes distance')+
  theme_Publication(base_family = 'sans')
print(g)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot.png', width = 5.5, height = 5.5, dpi = 300)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot.pdf', width = 5.5, height = 5.5, dpi = 300)
}

## random complex distance (use top 5% interaction proteins to caculate distance)
complex_distance_random_100 = list()
for(ni in 1:100){
  complex_distance_random = list()
  for(mi in 1:nrow(comp)){
    pro = strSplit(comp[mi, ], split = ' ')[1,]
    proNum = length(pro)
    proRandom = sample(proAll, proNum, replace=FALSE)
    compDist = di[proRandom, proRandom]
    upperDist = compDist[upper.tri(compDist)]
    dist_average = mean(upperDist)
    dist_min = min(upperDist)
    dist_max = max(upperDist)
    dist_median = median(upperDist)
    
    complex_distance_random[[mi]] = c(proNum = proNum, dist_average = dist_average, dist_min = dist_min, 
                                      dist_max = dist_max, dist_median = dist_median)
  }
  complex_distance_random_100[[ni]] = complex_distance_random
}
save(complex_distance_random_100, file = '../Rdata/19_00_complex_distance_random_100.Rdata')

## 
{
n = sample(1:100, 1, replace = FALSE)
dat = data.frame(do.call('rbind', complex_distance_random_100[[n]]))
datMelt = melt(dat, measure.vars = c('dist_average',   'dist_min',   'dist_max', 'dist_median'))
g = ggplot(data = datMelt, aes(variable, value)) +
  geom_boxplot() + 
  labs(x = '', y = 'Distance', title = 'Shuffled Complexes distance (one)')+
  theme_Publication(base_family = 'sans')
print(g)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_shuffled(one).png', width = 5.5, height = 5.5, dpi = 300)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_Shuffled(one).pdf', width = 5.5, height = 5.5, dpi = 300)
}

## caculate average of 100 random complexes sets
compRandom100 = lapply(complex_distance_random_100, function(x) {do.call('rbind', x)})

tmp = as.data.frame(do.call('rbind',compRandom100))
# tmp$random = rep(1:100, each = (nrow(tmp)/100))
tmp$comp = rep(paste0(rep('comp_', 5010), 1:5010), 100)
compRandom100_mean = lapply(unique(tmp$comp), function(mi){
  dat = subset(tmp, tmp$comp %in% mi)
  colAverage = colMeans(dat[,1:5])
  return(colAverage)
})
compRandom100_mean = do.call('rbind', compRandom100_mean)

## plot 
{
dat = data.frame(compRandom100_mean)
datMelt = melt(dat, measure.vars = c('dist_average',   'dist_min',   'dist_max', 'dist_median'))
g = ggplot(data = datMelt, aes(variable, value)) +
  geom_boxplot() + 
  labs(x = '', y = 'Distance', title = 'Shuffled Complexes distance (average of 100)')+
  theme_Publication(base_family = 'sans')
print(g)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_shuffled.png', width = 5.5, height = 5.5, dpi = 300)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_Shuffled.pdf', width = 5.5, height = 5.5, dpi = 300)
}


# dat = data.frame(compDist = complex_distance[,2], compdistRandom = complex_distance_random[,2])
## only take the mean dististance column
{
dat = data.frame(Predicted_Complexes_Distance = complex_distance[,2], 
                 Shuffled_Complexes_Distance = compRandom100_mean[,2])
datMelt = melt(dat)
g = ggplot(data = datMelt, aes(variable, value)) +
  geom_boxplot() + 
  labs(x = '', y = 'Distance', title = 'Complexes distance')+
  theme_Publication(base_family = 'sans', x_angle = -30)
# print(g)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_mean.png', width = 5, height = 7, dpi = 300)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_mean.pdf', width = 5, height = 7, dpi = 300)
}

### use one shuffled result for the mean distance column
{
n = sample(1:100, 1, replace = FALSE)
compRandom100 = data.frame(do.call('rbind', complex_distance_random_100[[n]]))

dat = data.frame(Predicted_Complexes_Distance = complex_distance[,2], 
                 Shuffled_Complexes_Distance = compRandom100[,2])
datMelt = melt(dat)
g = ggplot(data = datMelt, aes(variable, value)) +
  geom_boxplot() + 
  labs(x = '', y = 'Distance', title = 'Complexes distance (one)')+
  theme_Publication(base_family = 'sans', x_angle = -30)
# print(g)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_mean(one).png', width = 5, height = 7, dpi = 300)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot_mean(one).pdf', width = 5, height = 7, dpi = 300)
}
#

####### plot predicted complex and shuffled complex in one figure
## predicted complexes distance
predicted_complex = as.data.frame(complex_distance)
predicted_complexMelt = melt(predicted_complex, 
                             measure.vars = c('dist_average',  'dist_min',   'dist_max', 'dist_median'))
predicted_complexMelt$group = rep('Predicted', nrow(predicted_complexMelt))
## shuffeld complexes distance (mean of 100 times of shuffled complex sets)
dat = data.frame(compRandom100_mean)
datMelt = melt(dat, measure.vars = c('dist_average',   'dist_min',   'dist_max', 'dist_median'))
datMelt$group = rep('Shuffled', nrow(datMelt))

datPlot = rbind(predicted_complexMelt, datMelt)
g = ggplot(data = datPlot, aes(variable, value, fill = group)) +
  geom_boxplot() + 
  labs(x = '', y = 'Distance', title = 'Shuffled Complexes distance')+
  theme_Publication(base_family = 'sans')
print(g)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot(both).png', width = 6.5, height = 5.5, dpi = 300)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot(both).pdf', width = 6.5, height = 5.5, dpi = 300)


#################################
## random complex distance (use all proteins in MS/MS to caculate distance)
distance_all = get(load('../Rdata/19_00_distance_all.Rdata'))
proAll = rownames(distance_all)
complex_distanceAll_random_100 = lapply(1:100, function(ni){
  complex_distanceAll_random = lapply(1:nrow(comp), function(mi){
    pro = strSplit(comp[mi, ], split = ' ')[1,]
    proNum = length(pro)
    proRandom = sample(proAll, proNum, replace=FALSE)
    compDist = distance_all[proRandom, proRandom]
    upperDist = compDist[upper.tri(compDist)]
    dist_average = mean(upperDist)
    dist_min = min(upperDist)
    dist_max = max(upperDist)
    dist_median = median(upperDist)
    
    return(c(proNum = proNum, dist_average = dist_average, dist_min = dist_min, 
                                      dist_max = dist_max, dist_median = dist_median))
  })
   return(complex_distanceAll_random)
})
save(complex_distanceAll_random_100, file = '../Rdata/complex_distanceAll_random_100.Rdata')



## caculate average of 100 random complexes sets
compRandom100_all = lapply(complex_distanceAll_random_100, function(x) {do.call('rbind', x)})

tmp = as.data.frame(do.call('rbind',compRandom100_all))
# tmp$random = rep(1:100, each = (nrow(tmp)/100))
tmp$comp = rep(paste0(rep('comp_', 5010), 1:5010), 100)
compRandom100_all_mean = lapply(unique(tmp$comp), function(mi) {
  dat = subset(tmp, tmp$comp %in% mi)
  colAverage = colMeans(dat[,1:5])
  return(colAverage)
})
compRandom100_all_mean = do.call('rbind', compRandom100_all_mean)

####### plot predicted complex and shuffled complex in one figure
## predicted complexes distance
predicted_complex = as.data.frame(complex_distance)
predicted_complexMelt = melt(predicted_complex, 
                             measure.vars = c('dist_average',  'dist_min',   'dist_max', 'dist_median'))
predicted_complexMelt$group = rep('Predicted', nrow(predicted_complexMelt))
## shuffeld complexes distance (mean of 100 times of shuffled complex sets top 5% ppi proteins form MS/MS)
dat = data.frame(compRandom100_mean)
datMelt = melt(dat, measure.vars = c('dist_average',   'dist_min',   'dist_max', 'dist_median'))
datMelt$group = rep('Shuffled', nrow(datMelt))

## shuffeld complexes distance (mean of 100 times of shuffled complex sets all proteins form MS/MS)
dat2 = data.frame(compRandom100_all_mean)
datMelt2 = melt(dat2, measure.vars = c('dist_average',   'dist_min',   'dist_max', 'dist_median'))
datMelt2$group = rep('Shuffled_all', nrow(datMelt))


datPlot = rbind(predicted_complexMelt, datMelt, datMelt2)
g = ggplot(data = datPlot, aes(variable, value, fill = group)) +
  geom_boxplot() + 
  labs(x = '', y = 'Distance', title = 'Shuffled Complexes distance')+
  theme_Publication(base_family = 'sans')
print(g)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot(three).png', width = 6.5, height = 5.5, dpi = 300)
ggsave(g, filename = '../figs/19_00_complex_distance_boxplot(three).pdf', width = 6.5, height = 5.5, dpi = 300)


#########################

# distance plot for select complexes
# load('../Rdata/19_00_distance.Rdata')
distance = as.matrix(di)
# comp1
pro1 = comp[grep('MRPS12', comp[,1]),]
pro1 = t(strSplit(pro1[1], split = ' '))
pro1 = unique(pro1[,1])
# comp2 
pro2 = comp[grep('PSMA1', comp[,1]),]
pro2 = t(strSplit(pro2, split = ' '))
pro2 = unique(pro2[,1])
# comp3
pro3 = comp[grep('POLR3A', comp[,1]),]
pro3 = t(strSplit(pro3[1], split = ' '))
pro3 = unique(pro3[,1])
# comp4

pro4 = comp[grep('PML', comp[,1]),]
pro4 = t(strSplit(pro4[1], split = ' '))
pro4 = unique(pro4[,1])


# 
compPro = c(pro1, pro2, pro3, pro4)

# distance plot
compDist = distance[compPro, compPro]

# generate dataframe for annotation
annoRow = data.frame(Complex = c(rep('Complex1', length(pro1)),
                                 rep('Complex2', length(pro2)),
                                 rep('Complex3', length(pro3)),
                                 rep('complex4', length(pro4))))
rownames(annoRow) = compPro
p = pheatmap(compDist, color = colorRampPalette(c("firebrick3", "white", "navy"))(50), 
             border_color = NA, annotation_row = annoRow, cluster_cols = F, cluster_rows = F)
plotSave("../figs/19_00_CompDist_pheatmap.png", Plot = p, width = 7, height = 7)

# correlation plot
compAbun = compExp[compPro, ]
compCor = cor(t(compAbun), use = 'p', method = 'spearman')
p = pheatmap(compCor, color = colorRampPalette(c(  "navy", "white", "firebrick3"))(50), 
             border_color = NA, annotation_row = annoRow)
plotSave("../figs/19_00_CompCor_pheatmap.png", Plot = p, width = 7, height = 7)


# heatmap genes in order
dat = data.frame(index = 1:nrow(compDist), pro = rownames(compDist)[p$tree_row$order])