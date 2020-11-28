getSeqType <- function(line) {
  # projDir <- "~/Bishop.lab/Preprocessing/DRIP_Seq/GSE145964_Developmental_Context_ChIP_DRIP/"
  # setwd(projDir)
  # line <- "SRR11185367"
  jsonNow <- jsonlite::read_json(paste0("Data/QC/JSON/", line, ".json"))
  if("read2_before_filtering" %in% names(jsonNow)) {
    return("Yes")
  } else {
    return("No")
  }
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]

# Return result
result <- getSeqType(arg)
cat(result)