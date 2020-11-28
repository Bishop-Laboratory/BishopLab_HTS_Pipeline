parseConditions <- function(SRA, runInfoFile, groupFile, accList) {

  # SRA <- "A673_A5_STAG1"
  # projDir <- "~/Bishop.lab/Preprocessing/ChIP_Seq/Delattre_EWS_ChIPSeq_STAG2/"
  # setwd(projDir)
  # runInfoFile <- "Code/runInfo_2020.01.22_17.30.24.csv"
  # runInfo <- read.csv(runInfoFile, stringsAsFactors = F)
  # accList <- "Code/runInfo_2020.01.22_17.30.24.runList.txt"
  # groupFile <- "delattreGroupFile.txt"
  # 
  # groupInfo <- data.frame("SampleName" = unique(runInfo$SampleName),
  #                         "Group" = c(rep("CTR", 3), rep("2hr", 2), rep("24hr", 3)),
  #                         "Replicate" = 0,
  #                         "Condition" = c("S96", "Input", "RNH",
  #                                         "S96", "Input", "S96", "RNH", "Input"))
  # write.table(groupInfo, file = groupFile, sep = "\t", row.names = F, quote = F)
  
  # Proceed with merging of files
  
  # if(file.exists("Data/tmp/bamComparisonManifest_previousRun.RData")){
  #   load("Data/tmp/bamComparisonManifest_previousRun.RData")
  #   sraOldPrev <- sraOld
  #   
  # }
  # #availableBams <- list.dirs("/home/UTHSCSA/millerh1/DRIPSeqTest/Data/Bam_Files", recursive = F, full.names = F)
  # sraOld <- SRA
  # save(sraOld, file = "Data/tmp/bamComparisonManifest_previousRun.RData")
  
  availableBams <- list.dirs("Data/Bam_Files", recursive = F, full.names = F)
  runInfo <- read.csv(runInfoFile, fill = T, header = T, stringsAsFactors = F)
  groupInfo <- read.table(groupFile, sep = "\t", header = T, stringsAsFactors = F)
  if (! "PeakType" %in% colnames(groupInfo)) {
    groupInfo$PeakType <- "broad"
  }
  badSamps <- groupInfo[which(! groupInfo$SampleName %in% runInfo$SampleName),]
  
  if (! length(badSamps$SampleName)) {
    finalInfo <- merge(x = runInfo, y = groupInfo, by = "SampleName")
  } else {
    finalInfo1 <- merge(x = groupInfo, y = runInfo, by = "SampleName")
    finalInfo2 <- merge(x = badSamps, y = runInfo, by.y = "Experiment", by.x = "SampleName")
    finalInfo1 <- finalInfo1[,which(colnames(finalInfo1) %in% colnames(finalInfo2))]
    finalInfo2 <- finalInfo2[,which(colnames(finalInfo2) %in% colnames(finalInfo1))]
    finalInfo2 <- finalInfo2[,order(match(colnames(finalInfo2), colnames(finalInfo1)))]
    goodCols <- which(colnames(finalInfo1) == colnames(finalInfo2))
    finalInfo <- rbind(finalInfo1[,goodCols], finalInfo2[,goodCols])
  }
  if (length(finalInfo$SampleName) < 1) {
    return("Error: No samples available after merging.")
  }
  if ("SRAStudy" %in% colnames(finalInfo)) {
    if (all(finalInfo$Replicate == 0)) {
      finalInfo$resultName <- paste0(finalInfo$SRAStudy, "_", 
                                     finalInfo$Group, "_",
                                     finalInfo$Condition)
    } else {
      finalInfo$resultName <- paste0(finalInfo$SRAStudy, "_", 
                                     finalInfo$Group, "_",
                                     finalInfo$Condition, "_",
                                     finalInfo$Replicate)
    }
    
  } else {
    if (all(finalInfo$Replicate == 0)) {
      finalInfo$resultName <- paste0(
                                     finalInfo$Group, "_",
                                     finalInfo$Condition)
    } else {
      finalInfo$resultName <- paste0(
                                     finalInfo$Group, "_",
                                     finalInfo$Condition, "_",
                                     finalInfo$Replicate)
    }
  }
  accList <- read.table(accList, header = T)
  finalInfo <- finalInfo[order(match(finalInfo$Run, accList$Run)),] # Order by accession list
  if (! "Experiment" %in% colnames(finalInfo)) {
    finalInfo$Experiment <- finalInfo$SampleName
  }
  sampleNow <- finalInfo$SampleName[which(finalInfo$Run == SRA)]
  
  controlNow <- finalInfo$ControlSample[finalInfo$SampleName == sampleNow]
  alreadyRunGroups <- list.dirs("Results/bindingProfiles/", recursive = T, full.names = FALSE)
  alreadyRunGroups <- alreadyRunGroups[grep(alreadyRunGroups, pattern = paste0(availableBams, collapse = "|"))]
  prevInfo <- finalInfo[which(finalInfo$Run %in% availableBams & 
                                ! finalInfo$resultName %in% alreadyRunGroups),]
  if (is.na(controlNow)) {
    # CASE: You hit a control sample
    prevSamps <- prevInfo[which(prevInfo$ControlSample == sampleNow),]
    if (! length(prevSamps$SampleName)) {
      # CASE: You hit a control sample AND
      # no bam files exist for the experimental sample(s) that go with it 
      # OR a result has already been generated
      return("skip")
    } else {
      # CASE: You hit a control sample that does have bam files for it already in existence
      currentManifest <- prevSamps
      currentManifest$bamFile <- paste0("Data/Bam_Files/", currentManifest$Run, "/", currentManifest$Run, ".bam")
      controlRun <- SRA
      bamFileControl <- paste0("Data/Bam_Files/", controlRun, "/", controlRun, ".bam")
      currentManifest$comparisonString <- paste0(currentManifest$Group, "+",
                                                 currentManifest$resultName, "+", 
                                                 currentManifest$PeakType, "+",
                                                 currentManifest$bamFile, "+", bamFileControl)
      num <- sample(1:10000000, 1)
      file <- paste0("Data/tmp/bamComparisonManifest.list.", num, ".txt")
      currentManifest2 <- currentManifest[,which(colnames(currentManifest)== "comparisonString"), drop = FALSE]
      write.table(currentManifest2, quote = F, row.names = F, col.names = F,
                  file = file)
      return(file)
    }
  } else if (controlNow == "None") {
    # CASE: You hit an experimental sample for which no control exists
    currentManifest <- finalInfo[which(finalInfo$Run == SRA),]
    currentManifest$bamFile <- paste0("Data/Bam_Files/", currentManifest$Run, "/", currentManifest$Run, ".bam")
    controlRun <- "None"
    currentManifest$comparisonString <- paste0(currentManifest$Group, "+",
                                               currentManifest$resultName, "+", 
                                               currentManifest$PeakType, "+",
                                               currentManifest$bamFile, "+", controlRun)
    num <- sample(1:10000000, 1)
    file <- paste0("Data/tmp/bamComparisonManifest.list.", num, ".txt")
    currentManifest2 <- currentManifest[,which(colnames(currentManifest)== "comparisonString"), drop = FALSE]
    write.table(currentManifest2, quote = F, row.names = F, col.names = F,
                file = file)
    return(file)
  } else {
    # CASE: You hit an experimental sample
    prevSamps <- prevInfo[which(prevInfo$ControlSample == controlNow),]
    if (! length(prevSamps$SampleName)) {
      # CASE: You hit an experimental sample
      # BUT the corresponding control sample is not avaialble
      return("skip")
    } else {
      # CASE: You hit an experimental sample
      # AND the control sample is available
      currentManifest <- prevSamps
      currentManifest$bamFile <- paste0("Data/Bam_Files/", currentManifest$Run, "/", currentManifest$Run, ".bam")
      controlRun <- finalInfo$Run[finalInfo$SampleName == controlNow]
      bamFileControl <- paste0("Data/Bam_Files/", controlRun, "/", controlRun, ".bam")
      currentManifest$comparisonString <- paste0(currentManifest$Group, "+",
                                                 currentManifest$resultName, "+", 
                                                 currentManifest$PeakType, "+",
                                                 currentManifest$bamFile, "+", bamFileControl)
      num <- sample(1:10000000, 1)
      file <- paste0("Data/tmp/bamComparisonManifest.list.", num, ".txt")
      currentManifest2 <- currentManifest[,which(colnames(currentManifest)== "comparisonString"), drop = FALSE]
      write.table(currentManifest2, quote = F, row.names = F, col.names = F,
                  file = file)
      return(file)
    }
  }
  
  
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]
arg3 <- args[3]
arg4 <- args[4]

# Return result
result <- parseConditions(arg, arg2, arg3, arg4)
cat(result)
