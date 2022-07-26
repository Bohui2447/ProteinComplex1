

if(FALSE){
  library(funcTools)
  cleanComplex = function(compMCL){
    n = length(compMCL)
    stopifnot(is.list(compMCL))
    sublist = sapply(compMCL, is.list)
    if(!all(sublist)){
      if(!any(sublist)){
        return(compMCL[sapply(compMCL, length) > 1])
      }
      stop("No any element of compMCL is a list of complexes")
    }
    
    complex = list()
    for (mi in 1:n){
      complex = c(complex, compMCL[[mi]])
    }
    complex = unique(lapply(complex, sort))
    m = sapply(complex, length)
    return(complex[m > 1])
  }
  
  
  # predicted complex results
  if (FALSE){
    path = "../Rdata/17_00_MCL_onClusterOneResult/"
    file = dir(path, "*.Rdata")
    for(mi in file){
      fileName = paste0(path, mi)
      fileTxt = paste0("../Rdata/17_00_MCL_onClusterOneResult_txt/", tools::file_path_sans_ext(mi), ".txt")
      res = cleanComplex(get(load(fileName)))
      
      fid = file(fileTxt, "w")
      for(ni in res){
        cat(ni, file = fid)
        cat("\n", file = fid)
      }
      close(fid)
    }
  }
  
  # CORUM complex
  if (FALSE){
    path = "../data/"
    fileName = paste0(path, '10_CorumComp.Rdata')
    fileTxt = paste0("../Rdata/17_00_MCL_onClusterOneResult_txt/", tools::file_path_sans_ext(dir(path, "*.Rdata")), ".txt")
    res = get(load(fileName))
    res = res[, -1]
    fid = file(fileTxt, "w")
    for(ni in 1:nrow(res)){
      tmp = res[ni,]
      tmp = tmp[tmp != ""]
      cat(tmp, file = fid)
      cat("\n", file = fid)
    }
    close(fid)
  }
  
  path = "../Rdata/17_00_MCL_onClusterOneResult_txt/"
  pyFile = normalizePath("pyEvaluation/protein_complex_maps/complex_comparison.py")
  gold = normalizePath(paste0(path, '10_CorumComp.txt'))
  cmd = c()
  for(mi in dir(path, "*MCL\\.txt")){
    test = normalizePath(paste0(path, mi))
    cmd[mi] = paste0("cd /hpc/uu_npc/bli2/ProteinComplex/Rscript/pyEvaluation/protein_complex_maps; ", 
                     "python ", pyFile, " --cluster_predictions ", test, " --gold_standard ", gold)
  }
  
  ## qsub
  
  Rfile = thisFile()$path
  logFolder = paste0(thisFile()$dir, "/19_00_logFileFolder/")
  if(!dir.exists(logFolder)){
    dir.create(logFolder)
  }
  tmp = names(cmd)
  logFiles = paste0(logFolder, tmp, ".log")
  jobNames = paste0("complexEvalution_", tmp)
  path = getwd()
  rQsub2(path = path,
        rFile = Rfile,
        jobName = jobNames, threaded = 1,
        memoryG = 5, rTimeHour = 6, 
        email = "b.li1@uu.nl",
        when2Email = "FAIL,TIME_LIMIT",
        logFile = logFiles,
        preCMD = 'module load R/3.5.1 && Rscript ',
        param1 = cmd)
}

args <- commandArgs(TRUE)
system(args[1])
