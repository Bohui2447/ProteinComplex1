# 'Rdata_localPC' is 'RdataAll'

if(file.exists("../Rdata_localPC/4comProPro_CORUM.Rdata")){
  cat("load ../Rdata_localPC/4comProPro_CORUM.Rdata")
  load("../Rdata_localPC/4comProPro_CORUM.Rdata")
}else{
  source("00_00_function_complex.R")
  
  compHuman = read.csv(file = "../rawData_virtualBox/allComplexes_human.txt", 
                       header = T, sep = "\t", as.is = TRUE)
  # compHuman = read.csv(file = "../ComplexCORUM/allComplexes_human.txt", 
  #                      header = T, sep = "\t", as.is = TRUE)
  
  comp = data.frame(ComplexID = compHuman$ComplexID, 
                    geneName = as.character(compHuman$subunits.Gene.name.), 
                    stringsAsFactors = FALSE)
  proNum = sapply(comp$geneName, function(x) {
    y = strSplit(x, split = ";")
    n = ncol(y)
    return(n)})
  compNo1 = comp[proNum != 1,, drop = FALSE]
  colNm = c("ComplexID", paste0(paste0("V", 1:max(proNum)), collapse = ";"))
  # Using enough columns
  write.table(rbind(colNm, compNo1), file = "compNo1.csv",
              quote = FALSE, sep = ";", row.names = FALSE, col.names = FALSE)
  
  compNo1 = read.csv("compNo1.csv", sep = ";", header = TRUE)
  # save(compNo1, file = "../Rdata_localPC/10_CorumComp.Rdata")
  compNo1 = compNo1[rowSums(compNo1 != "") > 2,]
  
  # Complex-protein1-protein2
  comProPro = proComp(pm = compNo1, complexCol = 1)
  file.remove("compNo1.csv")
  # write.csv(comProPro, file = "../1rawData/comProPro.csv", row.names = FALSE)
  save(comProPro, file = "../Rdata_localPC/4comProPro_CORUM.Rdata")
}
