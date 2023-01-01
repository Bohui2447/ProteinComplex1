library(biomaRt)
library(stringr)
library(RCurl)
library(RJSONIO)
library(XML)
library(prideDownloader)
library(funcTools)


######## Uniprot ID to ENST
UniprotID2ENST = function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "uniprot_swissprot"
  # ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
  # ENSG00000163734       ENST00000296026    ENSP00000296026
  # attributes0 = c('uniprot_swissprot', 'entrezgene')
  attributes0 = c('uniprot_swissprot', 'ensembl_transcript_id')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}

######## ENSP to Uniprot ID 
ENSPToUniprotID= function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "ensembl_peptide_id"
  # ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
  # ENSG00000163734       ENST00000296026    ENSP00000296026
  # attributes0 = c('uniprot_swissprot', 'entrezgene')
  attributes0 = c('ensembl_peptide_id', "uniprot_swissprot")
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}


######## ENST to Uniprot ID 
ENST_To_UniprotID= function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "ensembl_transcript_id"
  # ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
  # ENSG00000163734       ENST00000296026    ENSP00000296026
  # attributes0 = c('uniprot_swissprot', 'entrezgene')
  attributes0 = c('ensembl_transcript_id', "uniprot_swissprot")
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}



######## Uniprot ID to gene name
UniprotID2geneName = function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "uniprot_swissprot"
  # ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
  # ENSG00000163734       ENST00000296026    ENSP00000296026
  # attributes0 = c('uniprot_swissprot', 'entrezgene')
  attributes0 = c('uniprot_swissprot', 'external_gene_name')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}

# atbt = listAttributes(mart) # get all mart in biomart
# write.csv(atbt, file = '../mart.csv')

#### extract useful data
extractData = function(PXD, colPatt = "Intensity.", proIDname = "Protein.IDs", namPXD, ref = "Gene"){
  # only extract Intensity.L. data
  # dataPXD_ = PXD[c(grep(proIDname, colnames(PXD)), grep(colPatt, colnames(PXD)))]
  
  # remove empty gene name, Contaminant and Reverse
  # namPXD$genes[is.na(namPXD$genes)] = ""
  # namPXD$Contaminant[is.na(namPXD$Contaminant)] = ""
  # namPXD$Reverse[is.na(namPXD$Reverse)] = ""
  # indxPXD = (namPXD$genes != "" & namPXD$Contaminant != "+" & namPXD$Reverse != "+")
  # dataPXD = data.frame(geneName = namPXD$genes[indxPXD], dataPXD_[indxPXD,])
  
  
  dataPXD = PXD
  geneNum = table(dataPXD[,ref])
  geneNam_larger1 = names(geneNum)[geneNum>1][-1]
  geneNam_equal1 = names(geneNum)[geneNum==1]
  
  # average duplicated rows
  dataPXD_noDupli0  = data.frame()
  t = 1
  len = length(geneNam_larger1)
  if(len > 0){
    for (x in geneNam_larger1){
      if (t %% 100 == 1) {cat(t, "/", len, "\n")}
      t = t+1
      x_data = dataPXD[dataPXD[,ref] == x,]
      if(ncol(x_data) < 4){
        x_data_mean = data.frame(as.character(unique(x_data[,ref])), 
                                 paste0(x_data[,ref], collapse = "|"), 
                                 mean(x_data[,3]), 
                                 stringsAsFactors = FALSE)
      }else{
        x_data_mean = data.frame(as.character(unique(x_data[,ref])), 
                                 paste0(x_data[,ref], collapse = "|"), 
                                 t(apply(x_data[, 3:ncol(x_data)], 2, function(x) mean(as.numeric(x), na.rm = T))), 
                                 stringsAsFactors = FALSE)
      }
      colnames(x_data_mean) = colnames(x_data)
      rownames(x_data_mean) = x
      dataPXD_noDupli0 = rbind(dataPXD_noDupli0, x_data_mean)
    }
  } else{
    dataPXD_noDupli0 = dataPXD
  }
  
  # combine duplicated and single rows
  dataPXD_equal1 = dataPXD[dataPXD[,ref] %in% geneNam_equal1, ]
  rownames(dataPXD_equal1) = dataPXD_equal1[,ref]
  dataPXD_noDupli = data.frame(rbind(dataPXD_noDupli0, dataPXD_equal1), 
                               stringsAsFactors = FALSE)
  
  return(dataPXD_noDupli)
}





##### combine PXD data into a big matrix
# dataPXDs = c("PXD_allData$PXD001553_catK_proteinGroups",  "PXD_allData$PXD001553_catL_proteinGroups" , "PXD_allData$PXD000194_proteinGroups")
combinePXDs = function(dataPXDs = c("dataPXD000309", "dataPXD000815", "dataPXD002619"), 
                       dataType = rep("Intensity.", 3), start = 13){
  n = 1
  colPXDall = geneNameAll = c()
  dataPXD0 = list()
  geneNameEach = list()
  for (mi in dataPXDs){
    dataPXD1 = eval(parse(text = mi))
    colPXD1 = sub(dataType[n], substr(mi, start, nchar(mi)), colnames(dataPXD1))
    colnames(dataPXD1) = colPXD1
    colPXDall = c(colPXDall, colPXD1[3:length(colPXD1)])
    # geneNameAll = c(geneNameAll, as.character(dataPXD1$Genes))
    geneNameEach[[n]] = as.character(dataPXD1$geneName)
    geneNameAll = c(geneNameAll, as.character(dataPXD1$geneName))
    dataPXD0[mi] = list(dataPXD1[, 3:length(colPXD1), drop = FALSE])
    n = n+1
  }
  geneNameAll = unique(geneNameAll)
  PXD_All = as.data.frame(matrix(0, nrow=length(geneNameAll), 
                                 ncol = length(colPXDall),
                                 dimnames = list(geneNameAll, colPXDall)))
  n = 1
  for(mi in dataPXDs){
    # PXD_All[rownames(dataPXD0[[mi]]), colnames(dataPXD0[[mi]])] = dataPXD0[[mi]]
    
    PXD_All[geneNameEach[[n]], colnames(dataPXD0[[mi]])] = dataPXD0[[mi]]
    
    n = n+1
  }
  return(PXD_All)
}





















################################################################## Downd PRIDE data Functions ####

################### count ###################
prideCount = function(query = "*",#query = "breast",
                      speciesFilter = NULL, #9606 is Homo sapiens
                      ptmsFilter = NULL,
                      tissueFilter = NULL,
                      diseaseFilter = NULL,
                      titleFilter = NULL,
                      instrumentFilter = NULL,
                      experimentTypeFilter = NULL,
                      quantificationFilter = NULL,
                      projectTagFilter = NULL){
  queryAll = unlist(as.list(environment()))
  cat(("your search query is: \n"))
  print(queryAll)
  cat("\n")
  queryAll = data.frame(type = names(queryAll), unname(queryAll))
  queryAll = paste(sapply(1:nrow(queryAll), function(x){
    paste(queryAll[x,1],queryAll[x,2], sep = "=")}), collapse  = "&")
  url = paste0('https://www.ebi.ac.uk/pride/ws/archive/project/count?', queryAll)
  # http://www.ebi.ac.uk/pride/ws/archive/project/count?query=breast&speciesFilter=9606ebi.ac.uk/pride/ws/archive/project/count?query=breast&speciesFilter=9606
  count = as.numeric(getURL(url))
  # count = as.numeric(suppressWarnings(readLines(url, 1)))
  cat("count = ", count, "\n")
  return(count)
}


################### search ###################
prideSearch = function(query = "*",#query = "breast",
                       speciesFilter = NULL, #9606 is Homo sapiens
                       ptmsFilter = NULL,
                       tissueFilter = NULL,
                       diseaseFilter = NULL,
                       titleFilter = NULL,
                       instrumentFilter = NULL,
                       experimentTypeFilter = NULL,
                       quantificationFilter = NULL,
                       projectTagFilter = NULL,
                       show = 100,#count # the maximum show is 2000
                       page = 0){
  queryAll = unlist(as.list(environment()))
  count = prideCount(query, speciesFilter, ptmsFilter, tissueFilter,
                     diseaseFilter,titleFilter, instrumentFilter,
                     experimentTypeFilter,quantificationFilter,
                     projectTagFilter)
  result1 = data.frame()
  for (i in 0:(ceiling(count/show) - 1)){
    queryAll["page"] = i
    queryAll = data.frame(type = names(queryAll), unname(queryAll))
    queryAll = paste(sapply(1:nrow(queryAll), function(x){
      paste(queryAll[x,1],queryAll[x,2], sep = "=")}), collapse  = "&")
    # url = paste0('http://www.ebi.ac.uk:80/pride/ws/archive/project/list?', queryAll)
    url = paste0('https://www.ebi.ac.uk/pride/ws/archive/project/list?', queryAll)
    # getting webpages
    text = getURL(url)
    if (length(grep("failure.html", text)) != 0){
      message("********** getURL fail **********")
      message("please try to reduce the number in query: 'show'!")
    } else{
      textJson = fromJSON(text)
      
      result = data.frame()
      for (mi in 1:length(textJson$list)){
        tmp = textJson$list[[mi]]
        result = rbind(result, data.frame(
          accession = tmp$accession,
          title = tmp$title,
          projectDescription = paste0(tmp$projectDescription, collapse = "|"),
          publicationDate = tmp$publicationDate,
          submissionType = tmp$submissionType,
          numAssays = tmp$numAssays,
          species = paste0(tmp$species, collapse = "|"),
          tissues = paste0(tmp$tissues, collapse = "|"),
          ptmNames = paste0(tmp$ptmNames, collapse = "|"),
          instrumentNames = paste0(tmp$instrumentNames, collapse = "|"),
          projectTags = paste0(tmp$projectTags, collapse = "|")))
      }
      result1 = rbind(result1, result)
    }
  }
  return(result1)
}

################### search for extra information of project ###################
prideSearchExtra = function(ProjectAccessionIDs){
  
  .searchExtra = function(ProjectAccessionID = "PXD008347"){
    url = paste0("https://www.ebi.ac.uk/pride/archive/projects/", ProjectAccessionID)
    html = getURL(url)
    doc = htmlParse(html, asText = TRUE)
    
    pathBase = "/html/body/div[2]/div/section/div/div[3]/div[2]/"
    apendix = c(Species = "div[1]/div[1]/p/a",
                Instrument = "div[3]/div[1]/p/a",
                Modification = "div[4]/div[1]/p/a",
                ExperimentType = "div[5]/div[1]/p/a",
                Tissue = "div[1]/div[2]/p/a",
                Software = "div[3]/div[2]/p/a",
                quantification = "div[4]/div[2]/p/a")
    
    xpath = paste0(pathBase, apendix)
    res = c()
    for(ni in seq_along(apendix)){
      tmp = paste0(xpathSApply(doc = doc, xpath[ni], xmlValue), collapse = "|")
      tmp = gsub("[\n][^A-Za-z]+", "", tmp)
      res = c(res, tmp)
    }
    res[res == ""] = "Not available"
    names(res) = names(apendix)
    return(res)
  }
  
  n = length(ProjectAccessionIDs)
  srchEx = matrix("", nrow = n, ncol = 7,
                  dimnames = list(ProjectAccessionIDs,
                                  names(.searchExtra(ProjectAccessionIDs[1]))))
  for(mi in seq_along(ProjectAccessionIDs)){
    m = ProjectAccessionIDs[mi]
    srchEx[m,] = .searchExtra(m)
    cat(mi, ":", m, "done\n")
    # Sys.sleep(1)
    if (((mi %% 10) == 0) && (mi != n)) {
      Sys.sleep(10)
      cat("PrideDatabase rule: wait 10 seconds")
    }
  }
  return(as.data.frame(srchEx))
}

#############################
xpath2ListID = function(xmlList,
                        xpath = "/html/body/div[2]/div/section/div/div[3]/div[2]/"){
  # ProjectAccessionID = "PXD006401"
  # url = paste0("https://www.ebi.ac.uk/pride/archive/projects/", ProjectAccessionID)
  # html = getURL(url)
  # doc = htmlParse(html, asText = TRUE)
  # xmlList = xmlToList(doc),
  
  l$body[[5]]$div$section$div[[6]][[6]][[2]][[2]]$p$a$text
  pathName = strsplit2(xpath, "/")[-(1:2)]
  
  # the number of each tag
  id = regexpr("\\d", pathName, perl=TRUE)
  id[id == -1] = 1
  tagID = regmatches(pathName, id)
  tagID[tagID == ""] = "1"
  tagID = as.numeric(tagID)
  tagID
  
  # the tag type
  tags = strsplit2(pathName, split = "\\[")[,1]
  
  # xpath to list ID
  tmp = xmlList
  n = NULL
  for(mi in seq_along(tags)){
    tmpName = names(tmp)
    tagMatch = which(tmpName == tags[mi])
    if (length(tagMatch) == 1){
      # print(tmpName[[tagMatch]])
      tmp = tmp[[tagMatch]]
      n = c(n, 1)
    } else {
      # print(tmpName[[tagMatch[tagID[mi]]]])
      tmp = tmp[[tagMatch[tagID[mi]]]]
      n = c(n, tagMatch[tagID[mi]])
    }
  }
  
  return(list(listID = n, res = tmp))
}

################### get links ###################
prideLink = function(ProjectAccessionID){
  result = data.frame()
  t = 1
  for (mi in ProjectAccessionID){
    systim = system.time({
      queryProjList = mi
      query = queryProjList
      # url = paste0('http://www.ebi.ac.uk:80/pride/ws/archive/file/list/project/', query)
      url = paste0('https://www.ebi.ac.uk/pride/ws/archive/file/list/project/', query)
      cat(paste(t,":", mi))
      textJson = fromJSON(getURL(url))
      cat(paste(" done, "))
      t = t+1
      for (ni in 1:length(textJson$list)){
        tmp = textJson$list[[ni]]
        result = rbind(result, data.frame(
          projectAccession = tmp$projectAccession,
          assayAccession = paste0(tmp$assayAccession, collapse = "|"),
          fileType = tmp$fileType,
          fileSource = tmp$fileSource,
          fileSize = tmp$fileSize,
          fileName = tmp$fileName,
          downloadLink = tmp$downloadLink))
      }
    })
    cat(paste("time:", unname(systim[3]), "\n"))
    if (t %% 8 == 0){
      cat("PrideDatabase rule: wait 5 seconds...\n")
      Sys.sleep(5)
    }
    if (t %% 29 == 0){
      cat("PrideDatabase rule: wait 10 seconds...\n")
      Sys.sleep(10)
    }
  }
  return(result)
}

prideDownload = function(prideLinks, fileMaxSize = 10, fileType = "SEARCH", path = getwd()){
  # prideLinks is the result from prideLink()
  link = prideLinks[(prideLinks$fileType %in% fileType),"downloadLink"]
  # what is the size (Gb) of these files
  totalSize = sum(prideLinks[(prideLinks$fileType %in% fileType), "fileSize"])/1000^3
  if (totalSize > fileMaxSize){
    stop("The size of total files is ",totalSize, " Gb > ", fileMaxSize,
         " Gb! \nPlease change fileMaxSize or reduce number of links")
  }
  cat("Total size of file is", totalSize,"Gb\n")
  if (substring(path, first = nchar(path)) != "/")
    path = paste0(path, "/")
  len = length(unique(link))
  t = 1
  for (ki in unique(link)){
    cat(t, "/",len,": Downloading", ki,"\n")
    t = t+1
    # tmp = strsplit(ki,"/")[[1]]
    # cmd = paste0('wget -q -O "', path,
    #              tmp[length(tmp)-1], "_", tmp[length(tmp)], '" "', ki, '"')
    # system(cmd)
    fileKi = paste0(path, basename(dirname(ki)), "_", basename(ki))
    download.file(ki, fileKi)
    if (t %% 10 == 0) Sys.sleep(abs(rnorm(1, mean = 10, sd = 5)))
  }
}

################################################################## Downd PRIDE data Functions ####







#################### functions copy from linux ###########

library(data.table)

proComp = function(pm, complexCol = 0, order = FALSE){
  if(complexCol == 0){
    compIDs = rownames(pm)
    proIDs = pm
  } else if (complexCol > 0){
    compIDs = pm[, complexCol, drop = TRUE]
    proIDs = pm[, -complexCol, drop = FALSE]
  } else {
    stop("complexCol must be non-negative !")
  }
  
  
  y = list()
  for(mi in 1:nrow(proIDs)){
    x = proIDs[mi,]
    x = x[x != ""]
    cb = t(combn(x, 2))
    y[[mi]] = cbind(compIDs[mi], cb)
  }
  
  comp = do.call(rbind, y)
  comp2 = comp[, c(1,3:2)]
  
  res = as.data.frame(rbind(comp, comp2), stringsAsFactors = F)
  
  res[, 2:3] = t(apply(res[, 2:3], 1, sort))
  res = res[!duplicated(res[, 2:3]),]
  if (order){
    res = as.data.frame(res[order(res[,2], res[,3]),])
  } else {
    res = as.data.frame(res)
  }
  rownames(res) = 1:nrow(res)
  
  colnames(res) = c("complex", "protein1", "protein2")
  return(res)
}



######## Uniprot ID to ENTREZ ID
library(biomaRt)
UniprotID2ENST = function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "uniprot_swissprot"
  # ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
  # ENSG00000163734       ENST00000296026    ENSP00000296026
  # attributes0 = c('uniprot_swissprot', 'entrezgene')
  attributes0 = c('uniprot_swissprot', 'ensembl_transcript_id')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}

## ensembl_gene_id to geneName
ENSGtoGeneName = function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "ensembl_gene_id"
  # ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
  # ENSG00000163734       ENST00000296026    ENSP00000296026
  # attributes0 = c('uniprot_swissprot', 'entrezgene')
  attributes0 = c('ensembl_gene_id', 'external_gene_name')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}

# atbt = listAttributes(mart) # get all mart in biomart
# write.csv(atbt, file = '../mart.csv')





## entrezgene to ensembl_gene_id
ENTREZtoENSG = function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "entrezgene"
  attributes0 = c('entrezgene', 'ensembl_gene_id')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}


## ensembl_gene_id to geneName
GeneNametoENSG = function(UniprotID, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "external_gene_name"
  # ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
  # ENSG00000163734       ENST00000296026    ENSP00000296026
  # attributes0 = c('uniprot_swissprot', 'entrezgene')
  attributes0 = c('external_gene_name', 'ensembl_gene_id')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}

##ensembl_gene_id to entrezgene_id
ENSGtoENTREZ = function(UniprotID, sameOrder = F){
  # host38.79 = "mar2015.archive.ensembl.org"
  # check www.ensembl.org/info/website/archives/index.html
  host38.89 = "may2017.archive.ensembl.org" 
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.89)
  # mart = useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host="www.ensembl.org")
  filters0 = "ensembl_gene_id"
  attributes0 = c(filters0,'entrezgene')
  # attributes0 = c(filters0,'entrezgene_id')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}

## entrezgene to geneName
ENTREZtoGeneName = function(UniprotID, sameOrder = F){
  host38.79 = "may2017.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = "entrezgene"
  attributes0 = c('entrezgene', 'external_gene_name')
  ID = getBM(filters = filters0, attributes = attributes0, values = UniprotID, mart = mart)
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[UniprotID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}




### Row sort # be careful: lowcase always bigger than uppercase letters
rowSort = function(mm){
  d <- data.table::as.data.table(mm)
  d[, row := .I]
  d <- melt(d, id.vars = "row") #wide to long format
  setkey(d, row, value) #sort
  d[, variable := paste0("V", ncol(mm):1)] #decreasing order
  
  #back to wide format and coerce to matrix
  msorted <- as.matrix(dcast(d, row ~ variable)[, row := NULL]) 
  return(msorted)
}


# F1 score
F1 = function(precision, recall){
  2*precision*recall/(precision+recall)
}

plot.prc = function(prc.perf, main = NULL, xlim = c(0, 1), ylim = c(0,1), 
                    plotF1 = TRUE, nudge_x = 0.02){
  # F1 = function(x) 2*x[1]*x[2]/(x[1]+x[2])
  
  data_prc = data.frame(prc.perf@x.values, prc.perf@y.values)
  colnames(data_prc) = c(prc.perf@x.name, prc.perf@y.name)
  if(is.na(data_prc[1,1])){
    data_prc[1,1] = 0
  } 
  if(is.na(data_prc[1,2])){
    data_prc[1,2] = 1
  } 
  # data_prc$F1 = apply(data_prc[,1:2], 1, F1)
  data_prc$F1 = F1(data_prc[,1], data_prc[,2])
  
  id = which.max(data_prc$F1)[1]
  data_F1 = data_prc[id,]
  g = ggplot(data_prc, aes_string(x = prc.perf@x.name, y = prc.perf@y.name))+
    geom_line()+
    geom_point(data = data_F1)+
    coord_cartesian(xlim = xlim, ylim = ylim)+
    ggtitle(main)+
    theme_Publication()
  if(plotF1) g = g +
    geom_text(aes(label = paste0("F1=", round(F1, 3))), data = data_F1, 
              hjust = "outward", nudge_x = nudge_x, color = "brown")
  return(g)
}





