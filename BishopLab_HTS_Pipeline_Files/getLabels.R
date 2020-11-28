getLabels <- function() {
  
  # bamFiles <- system('find ~/Bishop.lab/Preprocessing/DRIP_Seq/GSE81851_E2_MC/Data/Bam_Files/ -type f -name "*.bam"', 
  #                    intern = T)
  # sampleFile <- read.table("~/Bishop.lab/Preprocessing/DRIP_Seq/GSE81851_E2_MC/Code/downstreamSampleList.txt",
  #                          sep = "+", stringsAsFactors = F)
  
  sampleFile <- read.table("Code/downstreamSampleList.txt", sep = "+", stringsAsFactors = F)
  sampleFile$inputString <- gsub(pattern = "(^.*_.*)_.*$", 
                                 x = sampleFile$V2, replacement = "\\1_Input")
  bamFiles <- system('find Data/Bam_Files/ -type f -name "*.bam"', intern = T)
  # bamOrder <- gsub(pattern = ".*/Bam_Files/(.*)/.*$", x = bamFiles, replacement = "\\1")
  bamFrame <- data.frame(bamFiles = c(sampleFile$V3, sampleFile$V4),
                         fileLabels = c(sampleFile$V2, sampleFile$inputString))
  bamFrame <- unique(bamFrame)
  bamFrame <- bamFrame[order(match(bamFrame$bamFiles, bamFiles)),]
  labelsNow <- bamFrame$fileLabels
  output <- paste0(labelsNow, collapse = " ")
  return(output)
}

cat(getLabels())
