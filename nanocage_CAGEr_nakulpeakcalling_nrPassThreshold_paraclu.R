#---
#title: "Glioblastoma nanoCAGE Analysis"
#output: html_notebook
#---

#CAGEr Processing and expression analysis
#```{r}
library('Xmisc')
library("ggplot2")
library('CAGEr')
library('BSgenome.Hsapiens.UCSC.hg38')

args=commandArgs( trailingOnly=T )

idir=args[1]
ofile=args[2]

detach("package:Xmisc", unload=TRUE)

thread=16
gencodeReference <- "/bar/yliang/genomes/private/gencode.v36.annotation.gtf"

ifiles=list.files(path=idir, pattern="bam.nochrM.CTSS$", full.names=T)
print( length(ifiles) )

sampleLabels=gsub("-", "_", basename( ifiles )) #sample label cannot contain '-'. 
cs = new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg38", inputFiles = ifiles, inputFilesType = "ctss", sampleLabels = sampleLabels)

getCTSS(cs)
normalizeTagCount(cs, method="simpleTpm")

#corr.m <- plotCorrelation(cs, samples = "all", method = "pearson")
#nrPassThreshold=2 is from Brocks paper which may discard cell-line specific TSSs
clusterCTSS(object = cs, method="paraclu",  threshold = 0.1, nrPassThreshold = 1, thresholdIsTpm = TRUE, removeSingletons = TRUE, keepSingletonsAbove = 0.2, minStability = 2, maxLength = 100, reduceToNonoverlapping = TRUE, useMulticore = TRUE, nrCores = 4)

aggregateTagClusters(cs, tpmThreshold = 0.3, qLow = NULL, qUp = NULL, maxDist = 100, excludeSignalBelowThreshold = FALSE)

cc.df <- consensusClusters(cs, sample = NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL)
cc.df.tpm <- consensusClustersTpm(cs)

ctss.tags <- CTSStagCount(cs)

##
#Incorporation of unannotatedG peak composition
ifiles2=list.files(path=idir, pattern="unannotatedG.nochrM.CTSS$", full.names=T)
sampleLabels2=gsub("-", "_", basename( ifiles2 )) #sample label cannot contain '-'.
cs2 = new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg38", inputFiles = ifiles2, inputFilesType = "ctss", sampleLabels = sampleLabels2)

getCTSS(cs2)

ctss.tags.2 <- CTSStagCount(cs2)
##


library("snow")

clus <- makeCluster(thread)

clusterExport(clus, c("ctss.tags", "ctss.tags.2"))

resultapply <- parApply(clus,cc.df[,c("chr", "start", "end", "strand")],1,function(x) colSums(ctss.tags[ctss.tags$chr == x[1] & ctss.tags$pos >= as.numeric(x[2]) & ctss.tags$pos <= as.numeric(x[3]) & ctss.tags$strand == x[4], 4:ncol(ctss.tags)]))

resultapply2 <- parApply(clus,cc.df[,c("chr", "start", "end", "strand")],1,function(x) colSums(ctss.tags.2[ctss.tags.2$chr == x[1] & ctss.tags.2$pos >= as.numeric(x[2]) & ctss.tags.2$pos <= as.numeric(x[3]) & ctss.tags.2$strand == x[4], 4:ncol(ctss.tags.2)]))


stopCluster(clus)

#```

#Calculate the percentage unannotatedG reads for each peak. For now I use the max since if it is valid in one sample it should be valid in all

#```{r}
allCTSSpeakcount <- as.data.frame(t(resultapply))
uGCTSSpeakcount <- as.data.frame(t(resultapply2))

uGPercent <- uGCTSSpeakcount/allCTSSpeakcount

uGPercent[is.na(uGPercent)] <- 0

cc.df$maxUGPercentage <- apply(uGPercent, 1, function(x) max(x))

cc.df.tpm <- as.data.frame(cc.df.tpm)

cc.df <- cbind(cc.df, cc.df.tpm, uGPercent, allCTSSpeakcount, uGCTSSpeakcount)
#```

#Final filter 

#```{r}

#threshold <- .15

#cc.df.fil <- cc.df[cc.df$maxUGPercentage >= threshold,]

write.table(cc.df, ofile, quote=F, row.names=F, col.names=T )
#```
