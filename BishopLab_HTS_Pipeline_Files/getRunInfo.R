getRunInfo <- function(SRP, SRAList, downloadOnly, mode) {
  
  # # # # # # # Bug testing
  # SRP <- "none"
  # downloadOnly="None"
  # SRAList <- "Code/DRIP_accList.txt"
  # mode="ChIPSeq"
  # # SRAList <- "DRIP_accList_test.txt"

  toCheck <- c()
  if (SRP != 'none') {
    toCheck <- SRP
  }
  
  if (SRAList != 'none') {
    SRAList <- read.table(SRAList, header = F, stringsAsFactors = F)
    toCheck <- c(toCheck, SRAList$V1)
  }
  goodStr <- c("SRP", "SRR", "SRX", "SRS")
  toCheck2 <- toCheck
  toCheck <- toCheck[which(substr(toCheck, 1, 3) %in% goodStr)]
  if (! length(toCheck)) {
    SRAList2 <- SRAList
    colnames(SRAList2)[1] <- "Run"
    timevec <- Sys.time()
    timevec <- gsub(timevec, pattern =  ":", replacement =  ".")
    timevec <- gsub(timevec, pattern =  " ", replacement =  "_")
    timevec <- gsub(timevec, pattern =  "-", replacement =  ".")
    fileStr <- paste0("Code/runInfo_Table_", timevec, ".csv")
    #fileStr <- paste0("runInfo_Table_", timevec, ".csv")
    write.csv(SRAList2, file = fileStr, 
              quote = F, row.names = F, sep = ",")
    return(fileStr)
  }
  for (i in 1:length(toCheck)) {
    term <- toCheck[i]
    resp <- NULL
    while(is.null(resp)) {
      resp <- tryCatch(
        expr = {
          Sys.sleep(.5)
          httr::GET(url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi", 
                    query = list("save"="efetch","db"="sra","rettype"="runinfo", "term"=term))
        }, 
        error = function(cond) {
          Sys.sleep(.5)
          NULL
        }
      )
    }
    
    if (resp$status_code != 200) {
      msg <- paste0("ERROR: returned status code of ", 
                    as.character(resp$status_code),
                    " when contacting SRA API to obtain data for  ", term)
      return(msg)
    }
    if (httr::content(resp, as = "text") == "\n") {
      msg <- paste0("ERROR: Could not find valid entry for ", term)
      return(msg)
    }
    ct <- suppressMessages(read.csv(text = httr::content(resp, as = "text"), stringsAsFactors = F))
#    if (downloadOnly == "None") {
#      if (mode == "RNASeq") {
#        ct <- ct[which(ct$LibrarySource == "TRANSCRIPTOMIC"),]
#      } else if (mode == "ChIPSeq") {
#        ct <- ct[which(ct$LibrarySource == "GENOMIC"),]
#      }
#    }
    if (i == 1) {
      runInfo <- ct
    } else {
      cols <- which(colnames(runInfo) %in% colnames(ct))
      runInfo <- merge(x = runInfo, y = ct, by = cols, all = T)
    }
    
  }
  if (length(runInfo[,1]) == 0) {
    msg <- paste0("ERROR: No valid accessions provided -- no user input files detected either.",
                  " This error may also occur if names of pre-downloaded files start with 'SRR', 'SRX', or 'SRP' ",
                  " and are not valid SRA accessions.")
    return(msg)
  }
  # # Now get the group info and location info
  # groupInfo <- c()
  # locationInfo <- c()
  # queryVec <- runInfo$BioSample
  # for (i in 1:length(queryVec)) {
  #   term <- queryVec[i]
  #   resp <- NULL
  #   while(is.null(resp)) {
  #     resp <- tryCatch(
  #       expr = {
  #         Sys.sleep(.5)
  #         httr::GET(url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", 
  #                   query = list("db"="biosample", "id"=term))
  #       }, 
  #       error = function(cond) {
  #         NULL
  #       }
  #     )
  #   }
  #   if (resp$status_code != 200) {
  #     gi <- NA
  #     li <- NA
  #   } else {
  #     res <- XML::xmlToList(httr::content(resp, as = "text"))
  #     gi <- paste0(res$BioSample$Owner$Contacts$Contact$Name$First, " ",
  #                  res$BioSample$Owner$Contacts$Contact$Name$Last)
  #     li <- as.character(res$BioSample$Owner$Name)[1]
  #   }
  #   groupInfo <- c(groupInfo, gi)
  #   locationInfo <- c(locationInfo, li)
  # }
  n <- length(colnames(runInfo))-2
  runInfo2 <- runInfo[,c(1, 7, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 27, 28, 29:n)]
  runInfo2 <- as.data.frame(runInfo2)
  
  # Group info and location info are lists sometimes... Perhaps just merge with SRA info instead
  # groupInfo <- gsub(groupInfo, pattern = ",", replacement = " ")
  # locationInfo <- gsub(locationInfo, pattern = ",", replacement = " ")
  # 
  # runInfo2$group <- groupInfo
  # runInfo2$location <- locationInfo
  runInfo2 <- unique(runInfo2)
  timevec <- Sys.time()
  timevec <- gsub(timevec, pattern =  ":", replacement =  ".")
  timevec <- gsub(timevec, pattern =  " ", replacement =  "_")
  timevec <- gsub(timevec, pattern =  "-", replacement =  ".")
  fileStr <- paste0("Code/runInfo_", timevec, ".csv")
  #fileStr <- paste0("runInfo_Table_", timevec, ".csv")
  write.csv(runInfo2, file = fileStr, 
            quote = F, row.names = F)
  return(fileStr)
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]
arg3 <- args[3]
arg4 <- args[4]

# Return result
res <- suppressWarnings(suppressMessages(getRunInfo(arg, arg2, arg3, arg4)))
cat(res)

