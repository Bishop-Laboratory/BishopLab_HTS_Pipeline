#' DGE Pipeline
#'
#' Performs standard differential gene expression analysis using DESeq2
#'
#' @param pairedGenesList A list, named with primary genes of interest with
#'     vectors of secondary genes to test against OR a string containing the
#'     official MSIGDB name for a gene set of interest.
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#'     Either "All", "Normal_Tissues", or "Tumor_Tissues".
#'
#' @param outputPrefix Prefix for saved files. Should include directory info.
#'
#' @param sigTest Should the results be compared against random genes?
#'
#' @param nPerm Number of bootstrap sampling events to run during sigTest.
#'
#' @param plotMaxMinCorr If TRUE, the top correlated and anti-correlated genes
#'     will be plotted alongside the selected secondary genes.
#'
#' @param plotLabels If TRUE, correlation histograms will contain labeled lines showing
#'     secondary genes and their correlation values.
#'     If list of secondary genes is large, set this to FALSE or onlyTop to TRUE
#'     to avoid cluttering the plot.
#'
#' @param onlyTop For larger secondary gene lists -- This will filter the
#'     number of secondary genes which are plotted to avoid clutter if plotLabels = TRUE.
#'
#' @param topCutoff The value used for filtering if 'onlyTop' is 'TRUE'
#'
#' @param autoRug If the size of a secondary gene list > 50, plot lines will be replaced
#'     by an auto-generated rug. To disable this behavior, set to FALSE.
#'
#' @param returnDataOnly if TRUE will only return a list containing correlations
#'     and significance testing results if applicable.
#'
#' @return A list containing correlation values and signficance testing results
#'
#' @examples
#' pairedGenesList <- list("TP53" = c("BRCA1", "CDK12", "PARP1"),
#'                         "SON" = c("AURKB", "SFPQ", "DHX9"))
#'
#' pairedGenesList <- list("ATM" = "PUJANA_BRCA1_PCC_NETWORK")
#'
#' Result <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
#'                               Species = "hsapiens",
#'                               Sample_Type = "Normal_Tissues")
#'
#' @export
#'

DGE_pipeline <- function(
  projectDir,
  colData,
  contrast,
  toRun = c("Volcano", "circos", "GSEA", "Heatmap", "plotCounts", "dimReduce"),
  numTopHeat = 30,
  numTopPlotCount = 50,
  numCircosGenes = 2000,
  heatLegend = FALSE,
  dimReduceMethod = c("PCA", "distanceHeatmap"),
  returnDataOnly = FALSE,
  extraAnnoRows = NULL,
  dds = NULL,
  specDesign = NULL
) {

  # Bug testing
  projectDir <- "~/Bishop.lab/EWS_CTR_All_Cells/"
  colData <- "~/Bishop.lab/EWS_CTR_All_Cells/EWS_vs_MSC_colData.csv"


  # Preprocessing code



  # Libraries required
  require(EnhancedVolcano)
  require(clusterProfiler)
  require(pheatmap)
  require(tximport)
  require(circlize)
  require(data.table)
  require(ggpubr)
  require(biomaRt)
  require(RColorBrewer)
  require(gplots)
  require(DESeq2)

  if (! dir.exists(outdir) & ! returnDataOnly) {
    dir.create(outdir)
  }

  resList <- list()


  # Read in data
  contrastCol <- colnames(colData)[which(colnames(colData) == contrast[1])]
  # Wrangle colData
  colData <- colData[which(colData[, contrastCol] %in% c(contrast[2], contrast[3])),]
  colData[, contrastCol] <- factor(colData[, contrastCol], levels = c(contrast[3], contrast[2]))

  if (is.null(dds)) {
    print("Loading count data ... ")

    # Get count files
    samples <- rownames(colData)
    countFiles <- file.path(countDir, samples, "quant.sf")
    names(countFiles) <- samples
    # Load with Tximport
    txi <- tximport(files = countFiles,
                    type="salmon", tx2gene = tx2gene,
                    ignoreTxVersion = T)
    # Make DESeq dataset
    designStr <- paste0("~", contrastCol)
    if (! is.null(specDesign)) {
      designStr <- specDesign
    }
    dds <- DESeqDataSetFromTximport(txi = txi, colData = colData,
                                    design = formula(designStr))
    mn <- min(table(colData[, contrastCol]))
    keep <- rowSums(counts(dds)) >= mn
    dds <- dds[keep,]
    print("DESeq analysis ...")
    timestamp()
    dds <- DESeq(dds)
    print("DONE")
    timestamp()

    if (! returnDataOnly) {
      print("Saving dds object...")
      save(dds, file = file.path(outdir, "dds.RData"))
    }
    resList[["dds"]] <- dds
  } else {
    # Does DESeq need to be run still on user-supplied dds?
    gg <- resultsNames(dds)
    if (! length(gg)) {
      print("No results found in dds object -- running DESeq2...")
      timestamp()
      dds <- DESeq(dds)
      print("DONE")
      timestamp()
    }

    if (! returnDataOnly) {
      print("Saving dds object...")
      save(dds, file = file.path(outdir, "dds.RData"))
    }
    resList[["dds"]] <- dds

  }

  # Downstream of DESeq
  print("Gathering results from dds object ... ")
  #load("Results/DGE/EWS_vs_EWS/CHLA10_vs_CHLA9/dds.RData")
  # contrast <- c("Cell", "CHLA10", "CHLA9")
  res <- results(dds, contrast = contrast)
  # Change dds for plotting purposes
  dds <- dds[,which(colData(dds)[,contrastCol] %in% c(contrast[2], contrast[3]))]
  colData(dds)[,contrastCol] <- droplevels(colData(dds)[,contrastCol])

  resdf <- as.data.frame(res)
  resdf$geneName <- rownames(resdf)
  resdf <- resdf[,c(7,1:6)]
  resdf$padj[which(resdf$padj == 0)] <- .Machine$double.xmin
  resdf$GSEA <- -log10(resdf$padj) * sign(resdf$log2FoldChange)

  titleStr <- paste0(contrast[2], " vs. ", contrast[3])
  print("DONE")


  # Volcano
  if ("Volcano" %in% toRun) {
    print("Generating volcano plot ... ")
    pord <- resdf$padj
    pord <- -log10(pord)
    pord <- pord[order(pord, decreasing = T)]
    pval <- pord[150]
    ford <- resdf$log2FoldChange
    ford <- abs(ford)
    ford <- ford[order(ford, decreasing = T)]
    fval <- ford[600]
    ev <- EnhancedVolcano(resdf, lab = resdf$geneName, x = "log2FoldChange",
                          title = paste0(titleStr, " Volcano Plot"),
                          y = "padj", pCutoff = 10^-pval, FCcutoff = fval)
    if (! returnDataOnly) {
      ggsave(ev, filename = file.path(outdir, "Volcano.png"), height = 7, width = 8)
    }
    resList[["Volcano"]] <- ev
    print("DONE")
  }

  # Circos
  if ("circos" %in% toRun) {
    print("Generating circos plot ... ")
    cytopath <- file.path(tempdir(), 'ideoband.txt')
    cytopath2 <- paste0(cytopath, ".gz")
    if (species == "human") {
      cytobandURL <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz"
    }
    cmd <- paste0("wget -O ", cytopath2, " ", cytobandURL, " && gunzip -f ", cytopath2)
    system(cmd)
    cytoBand <- cytopath
    resdf2 <- resdf[, c(1, 8)]
    resdf2 <- unique(resdf2)
    pord <- abs(resdf2$GSEA)
    pord <- pord[order(pord, decreasing = T)]
    pval <- pord[numCircosGenes]

    res.df.sig <- resdf2[which(abs(resdf2$GSEA) >= pval),]

    if (is.null(mart.chr)) {
      ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
      mart.chr <- getBM(attributes = c("chromosome_name", "start_position", "end_position",
                                       "strand",
                                       "external_gene_name"),
                        mart = ensembl)
    }


    good.seq <- c(1:22, "X", "Y")
    mart.chr <- mart.chr[which(mart.chr$chromosome_name %in% good.seq),]
    mart.chr$chromosome_name <- paste0("chr", mart.chr$chromosome_name)
    colnames(mart.chr) <- c("seqname", "start", "end", "strand", "geneName")

    res.df.sig <- merge(x = mart.chr, y = res.df.sig, by = "geneName")
    res.df.sig2 <- res.df.sig[,c(2:4, 1, 6)]
    res.df.sig2 <- res.df.sig2[order(res.df.sig2$seqname, res.df.sig2$start),]
    qt <- quantile(abs(res.df.sig2$GSEA))
    topCol <- as.numeric(qt["75%"]) # use abs quartiles to determine color scale
    # topCol <- max(abs(res.df.sig2$GSEA)) * .85
    col_fun <- colorRamp2(c((-1 * topCol), 0, topCol), c("green", "black", "red"))

    resList[["circosData"]] <- res.df.sig2
    png(file.path(outdir, "DGE_circos.png"), height = 12, width = 12, units = c("in"), res = 600)
    circos.par(track.margin = c(.01, .01))
    circos.initializeWithIdeogram(cytoband = cytoBand, chromosome.index = paste0("chr", c(1:22, "X", "Y")))
    circos.genomicHeatmap(res.df.sig2, col = col_fun, numeric.column = 5)
    title(main = paste0("DGE circos ", titleStr))
    mtext(paste0("Top ", numCircosGenes, " significant genes"), side = 1)
    dev.off()
  }

  # Heatmap
  if ("Heatmap" %in% toRun) {
    print("Generating heatmap ... ")
    cts <- counts(dds, normalized = T)
    colData <- colData(dds)
    choose <- resdf[order(resdf$GSEA, decreasing = T),]
    choose <- choose[which(! is.na(choose$GSEA)),]
    n <- length(choose$geneName)
    l <- n-floor(numTopHeat/2)
    choose <- choose[c(1:floor(numTopHeat/2), l:n),]
    genes <- unique(as.character(choose$geneName))
    cts <- cts[genes,]
    colData <- as.data.frame(colData[which(as.character(colData[, contrast[1]]) %in% c(contrast[2], contrast[3])),])
    colDataPlt <- colData[, contrastCol, drop = F]
    if (! is.null(extraAnnoRows)) {
      colDataPlt <- colData[, c(contrastCol, extraAnnoRows), drop = F]
    }
    cts <- cts[,which(colnames(cts) %in% rownames(colDataPlt))]
    cts2 <- log2(cts + 1)
    cols <- c("lightcoral", "lightblue")
    names(cols) <- c(contrast[2], contrast[3])
    ph <- pheatmap(mat = cts2, annotation_col = colDataPlt,
                   scale = "row", color = greenred(100), legend = heatLegend,
                   show_colnames = F, main = paste0(titleStr, " DGE Heat Map\n") )
    if (! returnDataOnly) {
      ggsave(ph, filename = file.path(outdir, "Heatmap.png"), height = 8)
    }
    resList[["Heatmap"]] <- ph
    print("DONE")
  }

  # GSEA
  if ("GSEA" %in% toRun) {
    print("GSEA analysis ...")
    ranks <- resdf$GSEA
    names(ranks) <- resdf$geneName
    ranks <- ranks[which(! is.na(ranks))]
    ranks <- ranks[which(! duplicated(names(ranks)))]
    ranks <- ranks[order(ranks, decreasing = T)]
    GSEA <- myGSEA(ranks = ranks, returnDataOnly = returnDataOnly,
                   TERM2GENE = TERM2GENE,
                   outDir = outdir,
                   plotFile = "GSEA",
                   Condition = titleStr)
    resList[["GSEA"]] <- GSEA
    print("DONE")
  }

  # plotCounts
  if ("plotCounts" %in% toRun) {
    print("Plotting counts for top DGEs...")
    resdf2 <- resdf[order(resdf$GSEA, decreasing = T),]
    m <- floor(numTopPlotCount/2)
    n <- length(resdf2$geneName)
    geneToPlot <- resdf2[c(1:m, (n-m+1):n ) ,]
    k <- length(geneToPlot$geneName)
    plotCounts <- list()
    if (! returnDataOnly) {
      dir.create(file.path(outdir, "plotCounts"))
    }
    for (i in 1:k) {
      cat("\n", i, " of ", k)
      gene <- geneToPlot$geneName[i]
      subTitle <- paste0("FDR: ", signif(geneToPlot$padj[i], 3),
                         "\nLog2FoldChange: ", signif(geneToPlot$log2FoldChange[i], 3))
      # counts <- plotCounts(dds, gene = gene, intgroup = "condition", returnData = T)
      counts <- plotCounts(dds, gene = gene, intgroup = contrast[1], returnData = T)
      colnames(counts)[2] <- "condition"
      counts$count <- log2(counts$count + 1)
      gbp <- ggboxplot(data = counts, x = "condition", y = "count", add = "jitter",
                       ylab = "Normalized read counts", title = gene, subtitle = subTitle,
                       xlab = FALSE, fill = "condition") + rremove("legend")
      plotCounts[[i]] <- gbp
      names(plotCounts)[i] <- gene
      if (! returnDataOnly) {
        ggsave(gbp, filename = file.path(outdir, "plotCounts", paste0(gene, ".png")))
      }
    }
    resList[["plotCounts"]] <- plotCounts
  }

  # dimReduce
  if ("dimReduce" %in% toRun) {
    if ("PCA" %in% dimReduceMethod) {
      print("Calculating dimension reduction by PCA ... ")
      vsd <- vst(dds)
      resList[["vsd"]] <- vsd
      pca <- plotPCA(vsd, intgroup = contrast[1])
      resList[["pca"]] <- pca
      if (! returnDataOnly) {
        ggsave(pca, filename = file.path(outdir, "PCA.png"))
      }
      print("DONE")
    }
    if ("distanceHeatmap" %in% dimReduceMethod) {
      print("Calculating sample-to-sample distance from normalized read counts ... ")
      cts3 <- counts(dds, normalized = T)
      sampleDists <- dist(t(log2(cts3+1)))
      sampleDistMatrix <- as.matrix(sampleDists)
      colnames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
      anns <- as.data.frame(colData(dds))
      # annoCol <- anns[,"condition", drop = F]
      annoCol <- anns[,contrast[1], drop = F]
      if (! is.null(extraAnnoRows)) {
        print("Adding extra columns to heatmap")
        print(extraAnnoRows)

        annoCol <- anns[,c(contrast[1], extraAnnoRows), drop = F]
        print(annoCol)
      }
      ph <- pheatmap(sampleDistMatrix, annotation_row = annoCol,
                     clustering_distance_rows=sampleDists, show_rownames = F,
                     main = "Dim Reduction plot of gene expression\n",
                     clustering_distance_cols=sampleDists,
                     col=colors)

      resList[["dimReducePlot"]] <- ph
      resList[["dimReduceDistData"]] <- sampleDistMatrix
      resList[["dimReduceAnno"]] <- annoCol

      if (! returnDataOnly) {
        ggsave(ph, filename = file.path(outdir, "dimReducePlot.png"), height = 7, width = 7)
      }
      print("DONE")
    }


  }

  # Final output
  resdf <- resdf[which(! is.na(resdf$GSEA)),]
  cat("\nReturning results...\n")
  resdf <- merge(x = mart.res, y = resdf, by = "geneName")
  resdf <- resdf[order(resdf$padj),]
  if (! returnDataOnly) {
    fwrite(resdf, file = file.path(outdir, "DGE_results.csv"), sep = ",")
  }

  resList[["DGE_Result"]] <- resdf

  return(resList)

}

