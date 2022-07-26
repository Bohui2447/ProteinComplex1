

############# qsub 
if(FALSE){
  library(MCL)
  library(funcTools)
  
  options(max.print = 100)
  
  # num = length(param1)
  path = "/hpc/uu_npc/bli2/ProteinComplex/Rscript/"
  Rfile = thisFile()$path
  email = "b.li1@uu.nl"
  
  p1 = list.files('../data/', pattern = "*.tsv")
  p1 = factor(p1, unique(p1))
  p2 = as.character(c(1.2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15))
  
  
  params0 = expand.grid(p1, p2,  stringsAsFactors = T)
  params0$nppi = as.numeric(strSplit(strSplit(params0$Var1, split = "\\.")[,1], "_")[,4])
  params0 = sortDataframe(params0, by = "nppi", decreasing = TRUE)
  
  
  
  # id = paste0(apply(params0[,1:3], 1, paste0, collapse = "_"), "_ClusterOne.Rdata") %in% 
  #   dir("../Rdata/16_00_clusterOneResult_Rdata/")
  # 
  # 
  # 
  params = params0[ , 1:2]

  param1 = as.character(params[,1])
  param2 = as.character(params[,2])
  # param3 = as.character(params[,3])
  # param4 = as.character(params[,4])
  
  tmp = paste(param1, param2, sep = "_")
  logFolder = paste0(path, "18_00_logFileFolder/")
  if(!dir.exists(logFolder)){
    dir.create(logFolder)
  }
  logFiles = paste0(logFolder, tmp, ".log")
  jobNames = paste0("MCL_", tmp)
  
  # qsub Rjobs
  
  # res = rQsub(path = path, rFile = Rfile, jobName = jobNames[1:2],
  #             threaded = 1, rTimeHour = 10, memoryG = 20,
  #             logFile = logFiles[1:2], email = email,
  #             preCMD = "echo \"module load R/3.4.1 && Rscript ",
  #             param1 = param1[1:2], param2 = param2[1:2], 
  #             param3 = param3[1:2], param4 = param4[1:2])
  
  res = rQsub2(path = path, rFile = Rfile, jobName = jobNames,
               threaded = 1, rTimeHour = 30, memoryG = 30,
               logFile = logFiles, email = email, 
               when2Email = "FAIL,TIME_LIMIT",
               preCMD = "module load R/3.5.1 && Rscript ",
               param1 = param1, param2 = param2, 
               param3 = param3, param4 = param4)
  
  
  res
  View(qstat2())
  
}



path = "/hpc/uu_npc/bli2/ProteinComplex/data/"
args <- commandArgs(TRUE)
args1 = args[1]
ar1 = paste0(path, args[1]) # inputFile
ar2 = as.numeric(args[2]) # maxOverlap

# ar1 = paste0('../data/', param1[1])
# ar2 = as.numeric(param2[1])
# ar3 = as.numeric(param3[1])
# ar4 = as.numeric(param4[1])

print(.libPaths())
# print(installed.packages())
library('methods')
library('funcTools')
library('MCL')
# library('ClusterOneR')
print(ar1)
print(ar2)



##### load data 

# clusterResult = get(load(paste0(ar1, "_", ar2, "_", ar3, "_ClusterOne.Rdata")))
score = read.csv(paste0("/hpc/uu_npc/bli2/ProteinComplex/data/", args1),
                 header = T, as.is = T, sep = "\t")
# score = read.csv(paste0("/hpc/uu_npc/bli2/ProteinComplex/data/", params[1,1]),
#                  header = T, as.is = T, sep = "\t")

# generate score matrix
ProAll = unique(c(score$pro1, score$pro2))
scoreMat = edge2Mat(edge = score, rownm = ProAll, colnm = ProAll, direction = FALSE)
# generate complex matrix 

compMCL = mcl(x = compMat, addLoops = TRUE, 
          inflation = ar2, max.iter = 200)
compMCL = tapply(rownames(scoreMat), compMCL$Cluster, c)
  
folder = "/hpc/uu_npc/bli2/ProteinComplex/Rdata/18_00_MCL_onClusterOneResult/"
if(!dir.exists(folder)){
  dir.create(folder)
}


outPutFile = paste0(folder, args[1], "_", ar2, "_MCL.Rdata")

save(compMCL, file = outPutFile)






