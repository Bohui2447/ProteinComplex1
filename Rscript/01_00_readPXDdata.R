# 'Rdata_localPC' is 'RdataAll' 
# '../rawData_localPC/DownloadPride_20180827' is '../DownloadPride_20180827'

options(max.print = 200, stringsAsFactors = F)
rm(list = ls())
# setwd("E:/projects/ProteinComplexes_20180827/ProteinComplexes_script")
source("00_00_function_complex.R")

# read PXD***.txt and save PXD***.Rdata
file = dir(path = "../rawData_localPC/DownloadPride_20180827/LabelFree/", pattern = "*.txt")
fileAll = paste0("../rawData_localPC/DownloadPride_20180827/LabelFree/", file)
saveFile = paste0("../rawData_localPC/DownloadPride_20180827/LabelFree/LabelFreeRdata/", strSplit(file, "\\.txt")[,1], ".Rdata")

### The following samples are also included, because in the very initial analysis they were not excluded, 
### but actually they are not human samples. These should be excluded if further analysis will be done.
# toRM = c('PXD000661_organs_proteinGroups.txt', 'PXD000661_PrEST_proteinGroups.txt',
#          'PXD000661_tissues_proteinGroups.txt', 'PXD001548_C181_proteinGroups.txt',
#          'PXD001548_C182_proteinGroups.txt', 'PXD001548_C183_proteinGroups.txt',
#          'PXD001548_SAX1_proteinGroups.txt', 'PXD001548_SAX2_proteinGroups.txt',
#          'PXD001548_SAX3_proteinGroups.txt', 'PXD001548_SCX1_proteinGroups.txt',
#          'PXD001548_SCX2_proteinGroups.txt', 'PXD001548_SCX3_proteinGroups.txt',
#          'PXD001681_proteinGroups.txt', 'PXD002370_Ratio2_proteinGroups.txt',
#          'PXD002370_Ratio25_proteinGroups.txt', 'PXD002370_Ratio25Filter_proteinGroups.txt',
#          'PXD003422_proteinGroups.txt', 'PXD003531_PRIMARY_proteinGroups.txt',
#          'PXD003733_proteinGroups.txt', 'PXD003925_proteinGroups.txt',
#          'PXD003930_proteinGroups.txt', 'PXD005467_Host_proteinGroups.txt',
#          'PXD005467_TimecourseHostproteome_proteinGroups.txt', 'PXD005574_proteinGroups.txt',
#          'PXD006109_Cerebellum_proteinGroups.txt', 'PXD006109_DoE_proteinGroups.txt',
#          'PXD006109_MBR_proteinGroups.txt', 'PXD006109_noMBR_proteinGroups.txt',
#          'PXD006109_proteinGroups.txt', 'PXD006109_reproducibility_proteinGroups.txt',
#          'PXD006139_proteinGroups.txt', 'PXD006401_proteinGroups.txt',
#          'PXD006501_Main_proteinGroups.txt', 'PXD006501_Opt_proteinGroups.txt',
#          'PXD006579_proteinGroups.txt', 'PXD006897_proteinGroups.txt',
#          'PXD007039_proteinGroups.txt', 'PXD008289_proteinGroups.txt',
#          'PXD008691_HEK_proteinGroups.txt', 'PXD008691_neuron_proteinGroups.txt',
#          'PXD008720_DDA_proteinGroups.txt', 'PXD009792_proteinGroups.txt')
# toRM1 = paste0("../rawData_localPC/DownloadPride_20180827/LabelFree/", toRM)
# fromRM = paste0("../rawData_localPC/DownloadPride_20180827/LabelFree/removed/", toRM)
# if (!all(file.exists(toRM1))){
#   file.copy(fromRM, toRM1)
# }
# file.remove(toRM1)

if(!all(file.exists(saveFile))){
  for (mi in (1:length(saveFile))){
    if (file.exists(saveFile[mi])){
      load(saveFile[mi])
    } else {
      print(mi)
      print(fileAll[mi])
      
      if(fileAll[mi] == "../rawData_localPC/DownloadPride_20180827/LabelFree/PXD000608_proteinGroups.txt"){
        quote = "\""
      } else {
        quote = ""
      }
      PXD = read.csv(file = fileAll[mi], sep = "\t", header = T, quote = quote, 
                     na.strings=c("NA", "NaN"), stringsAsFactors=F, comment.char = "")
      save(PXD, file = saveFile[mi])
    }
  }
}

# save 1PXD_all.Rdata
PXD_all = list()
fname = strSplit(strSplit(saveFile, split = "/")[,6], split = "\\.")[,1]
for(mi in seq_along(saveFile)){
  PXD_all[[fname[mi]]] = get(load(saveFile[mi]))
}
save(PXD_all, file = "../Rdata_localPC/0PXD_all.Rdata")



