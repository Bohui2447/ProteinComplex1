# 'Rdata_localPC' is 'RdataAll'


rm(list = ls())
# setwd("E:/projects/ProteinComplexes_20180808/ProteinComplexes_script")
source("00_00_function_complex.R")
PXD_intensityGeneName = get(load("../Rdata_localPC/2PXD_intensityGeneName.Rdata"))

#### remove geneName with NA
PXD_intensityGeneName = lapply(PXD_intensityGeneName, function(x){
  # y = x[,-2] # remove protein ID
  y = x[!is.na(x[,1]), ]
  # rownames(y) = y[,1]
  return(y)
})


#### combine PXD_data
if (file.exists("../Rdata_localPC/3PXD_combine.Rdata")){
  cat("load 3PXD_combine.Rdata")
  x = load("../Rdata_localPC/3PXD_combine.Rdata")
} else {
  PXD_allData = list()
  for (mi in names(PXD_intensityGeneName)){
    print(mi)
    PXD = PXD_intensityGeneName[[mi]]
    namPXD = PXD[, 1:2]
    namPXD = namPXD[, c(2,1)]
    colnames(namPXD) = c("uniprotIDs", "geneName")
    # need to change ***** colPatt and Protein.IDs ************
    PXD_allData[[mi]] = extractData(PXD = PXD, colPatt = "Intensity.", 
                                    proIDname = "Protein.IDs", namPXD = namPXD, ref = "geneName")
  }
  PXD_combine = combinePXDs(dataPXDs = paste("PXD_allData", names(PXD_allData), sep = "$"), 
                            dataType = rep("Intensity", length(PXD_allData)), start = 13)
  sum(is.na(PXD_combine))
  save(PXD_combine, file = "../Rdata_localPC/3PXD_combine.Rdata")
}

View(t(t(sort(colSums(PXD_combine)))))
View(t(PXD_combine[1:5,]))



##############################


# # How does number of proteins increase by datasets
# row = c()
# for (ni in 1:74){
#   PXD_combine = combinePXDs(dataPXDs = paste("PXD_allData", names(PXD_allData)[1:ni], sep = "$"), 
#                             dataType = rep("Intensity", ni), start = 13)
#   row = c(row, nrow(PXD_combine))
# }
# png("../How does number of proteins increase by datasets.png", 
#     units = "in", width = 6, height = 6, res = 300)
# plot(1:74, row, ylim = c(0, 1.1*max(row)), xlab = "number of datasets", 
#      ylab = "number of proteins", main = "How does number of proteins increase by datasets" )
# dev.off()



