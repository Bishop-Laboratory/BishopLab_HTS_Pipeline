get_SYN_info <- function(synID, synList, scriptDIR) {
  
  # # # # # # Bug testing
  # synID <- "none"
  # synList <- "/home/UTHSCSA/millerh1/Bishop.lab/Preprocessing/RNA_Seq/ROSMAP_DLPFC/finalAccList.txt"
  # scriptDIR <- "~/bin/BishopLab_HTS_Pipeline_Files/"
  dbdat <- read.csv(file.path(scriptDIR, "synapse_studies_list.csv"), stringsAsFactors = F)
  dbdatFinal <- NULL
  if (synList != 'none') {
    synList <- read.table(synList, header = F, stringsAsFactors = F)
    
    dbdatFinal <- dbdat[which(dbdat$specimenID %in% synList[,1] | 
                            dbdat$name %in% synList[,1] | 
                            dbdat$id %in% synList[,1]),]
    
  } 
  if (synID != "none") {
    dbdat3 <- dbdat[which(dbdat$parentId == synID | dbdat$projectId == synID),]
    if (! is.null(dbdatFinal)) {
      dbdatFinal <- dbdat3
    } else {
      dbdatFinal <- rbind(dbdatFinal, dbdat3)
    }
  }
  
  
  dbdatFinal <- unique(dbdatFinal)
  colnames(dbdatFinal)[3] <- "Run"
  timevec <- Sys.time()
  timevec <- gsub(timevec, pattern =  ":", replacement =  ".")
  timevec <- gsub(timevec, pattern =  " ", replacement =  "_")
  fileStr <- paste0("Code/runInfo_", timevec, ".csv")
  write.csv(dbdatFinal, file = fileStr, 
              quote = F, row.names = F)
  return(fileStr)
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]
arg3 <- args[3]

# Return result
res <- suppressWarnings(suppressMessages(get_SYN_info(arg, arg2, arg3)))
cat(res)

