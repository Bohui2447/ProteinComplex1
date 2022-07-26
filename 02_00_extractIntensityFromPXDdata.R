# 'Rdata_localPC' is 'RdataAll'
# '../rawData_localPC/DownloadPride_20180827' is '../DownloadPride_20180827'

rm(list = ls())
source("00_00_function_complex.R")
load("../Rdata_localPC/0PXD_all.Rdata")
# remove proteinID is IPI stat
# index = c("PXD000254_proteinGroups","PXD001487_proteinGroups","PXD003095_proteinGroups",
#           "PXD003438_HuH7_proteinGroups","PXD003438_N2a_proteinGroups","PXD005559_proteinGroups",
#           "PXD005940_proteinGroups","PXD005942_proteinGroups","PXD005946_proteinGroups",
#           "PXD007725_proteinGroups", "PXD007203_proteinGroups", "PXD006722_proteinGroups",
#           "PXD006775_proteinGroups", "PXD001430_4AB_proteinGroups", "PXD004980_5305_proteinGroups",
#           "PXD004980_5646_proteinGroups", "PXD004980_5646db_proteinGroups","PXD004980_5865_proteinGroups")
# for(ni in index){
#   PXD_all[ni] = NULL
# }


#### extract data from PXD_all
if (file.exists("../Rdata_localPC/1PXD_intensity.Rdata")){
  cat("load ../Rdata_localPC/1PXD_intensity.Rdata\n")
  load("../Rdata_localPC/1PXD_intensity.Rdata")
} else{
  datasetInfo = read.csv("../rawData_localPC/DownloadPride_20180827/LabelFree/datasetInfo.csv", 
                         row.names = 1, stringsAsFactors = FALSE) # 
  
  PXD_intensity0 = list()
  for(mi in names(PXD_all)){
    print(mi)
    PXD = PXD_all[[mi]]
    # extract Protein.IDs
    # Protein.IDs = as.character(PXD[, grep(datasetInfo[mi, "proIDname"], x = colnames(PXD), 
    #                                       value = TRUE, ignore.case = FALSE, fixed = T)])
    
    # if (!mi %in% rownames(datasetInfo)){
    #   next
    # }
    
    Protein.IDs = as.character(PXD[, datasetInfo[mi, "proIDname"]])
    
    if (datasetInfo[mi, "process"] == "split Protein.IDs"){
      # tmp = str_extract_all(Protein.IDs, pattern = "\\|\\w*\\|")
      tmp = str_extract_all(Protein.IDs, pattern = "\\|\\w*\\|")
      Protein.IDs = sapply(tmp, function(x) {paste0(sapply(x, function(y) {
        substr(y, start = 2, stop = nchar(y)-1)}),
        collapse = ";")})
    }
    
    # deal with PXD003935_MP1to12_proteinGroups 
    if(mi == "PXD003935_MP1to12_proteinGroups"){
      Protein.IDs = strSplit(Protein.IDs, split = "\\|")[,1]
    }
    
    # deal with PXD007740_proteinGroups 
    if(mi == "PXD007740_proteinGroups"){
      Protein.IDs = strSplit(Protein.IDs, split = "\\_")[,1]
    }
    
    
    # deal with ensP proteinIDs
    ind = c("PXD000895_proteinGroups",
            "PXD006138_proteinGroups", "PXD008047_proteinGroups")
    if(mi %in% ind){
      ensP = strSplit(Protein.IDs, split = "_")[,1]
      tmp = ENSPToUniprotID(ensP)
      proID = sapply(ensP, function(x) {
        if(x %in% tmp$ensembl_peptide_id){
          uniprot = tmp[tmp$ensembl_peptide_id == x, 2][1]
        } else 
          uniprot = ""
        return(uniprot)})
      Protein.IDs = unname(unlist(proID))
    }
    
    # deal with ensT proteinIDs
    if(mi == "PXD005445_proteinGroups"){
      ensP = strSplit(Protein.IDs, split = "\\.")[,1]
      tmp = ENST_To_UniprotID(ensP)
      proID = sapply(ensP, function(x) {
        if(x %in% tmp$ensembl_transcript_id){
          uniprot = tmp[tmp$ensembl_transcript_id == x, 2][1]
        } else 
          uniprot = ""
        return(uniprot)})
      Protein.IDs = unname(unlist(proID))
    }
    
    # deal with UniRef100 and CONTAMINANT_sp, "CON__"
    if(mi %in% c("PXD003134_mHyo_proteinGroups", "PXD003133_100min_proteinGroups",
                 "PXD003134_proteinGroups", "PXD003133_20min_proteinGroups",
                 "PXD007686_proteinGroups")){
      Protein.IDs[grep("CON__", Protein.IDs)] = ""
      Protein.IDs[grep("UniRef100", Protein.IDs)] = ""
      Protein.IDs[grep("CONTAMINANT_sp", Protein.IDs)] = ""
    }
    
    ## deal with dataset tr|***
    specInd = c("PXD006486_IsoSeq_proteinGroups", "PXD007686_proteinGroups",
                "PXD006486_siSF3B1_proteinGroups", "PXD006486_siSRSF1_proteinGroups",
                "PXD008260_proteinGroups")
    if (mi %in% specInd){
      # tmp = str_extract_all(Protein.IDs, pattern = "\\|\\w*\\|")
      tmp = str_extract_all(Protein.IDs, pattern = "sp\\|\\w*\\|")
      Protein.IDs = sapply(tmp, function(x) {paste0(sapply(x, function(y) {
        substr(y, start = 4, stop = nchar(y)-1)}),
        collapse = ";")})
    }
    
    Protein.IDs = strSplit(Protein.IDs, split = ";")[,1]
    # Protein.IDs have ###CONTAMINANT###
    specialID = grep("###", Protein.IDs)
    if(length(specialID) > 0){
      Protein.IDs[specialID] = ""
    }
    
    # extract Reverse
    Reverse = PXD[, grep(datasetInfo[mi, "Reverse"], x = colnames(PXD), 
                         value = TRUE, ignore.case = TRUE)]
    Reverse[is.na(Reverse)] = ""
    # extract Contaminant
    Contaminant = PXD[, grep(datasetInfo[mi, "Contaminant"], x = colnames(PXD), 
                             value = TRUE, ignore.case = TRUE)]
    Contaminant[is.na(Contaminant)] = ""
    # extract peptideCount
    PeptideCount = PXD[, grep(datasetInfo[mi, "PeptideCount"], x = colnames(PXD), 
                              value = TRUE, ignore.case = TRUE)]
    PeptideCount = strSplit(PeptideCount, split = ";")[,1]
    PeptideCount = as.numeric(PeptideCount)
    PeptideCount[is.na(PeptideCount)] = 0
    # extract intensity
    intensity = PXD[, grep(datasetInfo[mi, "LightOrHeavy"], x = colnames(PXD), 
                           value = TRUE, ignore.case = TRUE), drop = FALSE]
    PXD_inten = data.frame(Protein.IDs, PeptideCount, Reverse, Contaminant, 
                           intensity, stringsAsFactors = F)
    PXD_intensity0[[mi]] = PXD_inten
  }
  
  
  
  # remove low quality data
  PXD_intensity = lapply(PXD_intensity0, function(x){
    y = subset(x, (Contaminant != "+") & (Reverse != "+") & (PeptideCount > 1))
    y[,-(2:4)]
  })
  save(PXD_intensity, file = "../Rdata_localPC/1PXD_intensity.Rdata") 
}


# protein_ID to gene name
if(file.exists("../Rdata_localPC/2PXD_intensityGeneName.Rdata")){
  cat("load ../Rdata_localPC/2PXD_intensityGeneName.Rdata\n")
  load("../Rdata_localPC/2PXD_intensityGeneName.Rdata")
} else {
  
  proID = unique(unlist(lapply(PXD_intensity, function(x) x$Protein.IDs)))
  proID = proID[proID != ""]
  uniprotGeneName = UniprotID2geneName(UniprotID = proID)
  
  uniprotGeneName = uniprotGeneName[!duplicated(uniprotGeneName$uniprot_swissprot),]
  rownames(uniprotGeneName) = uniprotGeneName$uniprot_swissprot
  PXD_intensityGeneName = lapply(setNames(names(PXD_intensity), names(PXD_intensity)), function(mi){
    x = PXD_intensity[[mi]]
    y = x[x$Protein.IDs != "",]
    if (mi == "PXD001170_PTP_proteinGroups"){
      y$geneName = y$Protein.IDs
    } else {
      y$geneName = uniprotGeneName[y$Protein.IDs,2]
    }
    return(y[, c(ncol(y), 1:(ncol(y)-1))])
  })
  # PXD_intensityGeneName = lapply(PXD_intensity, function(x){
  #   y = x[x$Protein.IDs != "",]
  #   y$geneName = uniprotGeneName[y$Protein.IDs,2]
  #   return(y[, c(ncol(y), 1:(ncol(y)-1))])
  # })
  
  save(PXD_intensityGeneName, file = "../Rdata_localPC/2PXD_intensityGeneName.Rdata")
}