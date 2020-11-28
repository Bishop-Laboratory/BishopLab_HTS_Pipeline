techReps <- function(line, runInfoFile) {
  # projDir <- "~/Bishop.lab/Preprocessing/DRIP_Seq/RLoopDB_All_Available_DRIP_Data/"
  # setwd(projDir)
  # runInfoFile <- "Code/runInfo_2020.01.08_18.43.31.csv"
  # line <- "SRR10168931"

  # Proceed with merging of files
  runInfo <- read.csv(runInfoFile, fill = T, header = T, stringsAsFactors = F)
  if (line %in% runInfo$Run) {
    exp <- runInfo$Experiment[which(runInfo$Run == line)]
    expDups <- runInfo$Experiment[which(duplicated(runInfo$Experiment))]
    if (exp %in% expDups) {
      res <- 'yes'
      lines <- runInfo$Run[which(runInfo$Experiment == exp)]
      ind <- which(lines == line)
      if (ind == 1) {
        lines <- as.matrix(lines)
        write.table(lines, file = "Data/tmp/runTempTable.txt", quote = F, row.names = F)
      } else {
        res <- "skip"
      }
    } else {
      res <- 'no'
    }
  } else {
    res <- 'no'
  }
  return(res)
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]

# Return result
result <- techReps(arg, arg2)
cat(result)