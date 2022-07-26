# run ClusterOne on ppi_predict (different score)

if(FALSE){
  library(funcTools)
  library(ClusterOneR)
  Rfile = thisFile()$path
  {
    # setwd2thisFile()
    options(max.print = 100)
    
    path = "/hpc/uu_npc/bli2/ProteinComplex/Rscript"
    
    email = "b.li1@uu.nl"
    
    p1 = list.files('../data/', pattern = "*.tsv")
    p1 = factor(p1, unique(p1))
    p2 = as.character(c(0.2, 0.25, 0.3, 0.35, 0.4))
    p3 = c(0.6, 0.7, 0.8)
    
    params0 = expand.grid(p1, p2, p3, stringsAsFactors = T)
    params0$nppi = as.numeric(strSplit(strSplit(params0$Var1, split = "\\.")[,1], "_")[,4])
    params0 = sortDataframe(params0, by = "nppi", decreasing = TRUE)
    
    id = paste0("clusterResult_", apply(params0, 1, paste0, collapse = "_"), 
                ".Rdata") %in% dir("../Rdata/16_00_clusterOneResult_Rdata/")
    params = params0[!id, 1:3]
    
    param1 = as.character(params[,1])
    param2 = as.character(params[,2])
    param3 = as.character(params[,3])
    
    tmp = apply(params, 1, function(x)paste0(x, collapse = "_"))
    
    logFolder = paste0(path, "/16_00_logFileFolder/")
    if(!dir.exists(logFolder)){
      dir.create(logFolder)
    }
    logFiles = paste0(logFolder, tmp, ".log")
    jobNames = paste0("networkClusterOne_", tmp)
    
    res = rQsub2(path = path, rFile = Rfile, jobName = jobNames,
                 threaded = 1, rTimeHour = 2, memoryG = 30,
                 logFile = logFiles, email = email,
                 preCMD = "module load R/3.5.1 && Rscript ",
                 param1 = param1, param2 = param2, param3 = param3)
    
    res
    View(qstat2())
  }
}

path = "/hpc/uu_npc/bli2/ProteinComplex/data/"
args <- commandArgs(TRUE)
ar1 = paste0(path, args[1]) # inputFile
ar2 = as.numeric(args[2]) # minDensity
ar3 = as.numeric(args[3]) # maxOverlap
print(.libPaths())
# print(installed.packages())
library('methods')
library('funcTools')
library('ClusterOneR')
print(ar1)
print(ar2)
print(ar3)

clusterResult = clusterOneR(inputFile = ar1, minDensity = ar2, minSize = 2,
                            maxOverlap = ar3, seedMethod = "nodes")
folder = "/hpc/uu_npc/bli2/ProteinComplex/Rdata/16_00_clusterOneResult_Rdata/"
if(!dir.exists(folder)){
  dir.create(folder)
}
outPutFile =  paste0(folder, args[1], "_", ar2, "_", ar3, "_ClusterOne.Rdata")

save(clusterResult, file = outPutFile)

print(mem_used2(T))



