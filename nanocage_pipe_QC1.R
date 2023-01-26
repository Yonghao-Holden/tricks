#!/usr/bin/env Rscript
#####################################################################
# This is a Rscript to plot reads and peaks distribtion of nanoCAGE data
#####################################################################
library('Xmisc')
library("data.table")
library("preprocessCore")
library("dplyr")

# Default parameters

# User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-wd',type='character',help='working directory')
parser$add_argument('-input',type='character',help='input filtered bam file')
parser$add_argument('-samplename',type='character',help='sample name')
parser$add_argument('-genome',type='character',help='genome')
argsparse <- parser$get_args()

# Argument parse to change defaults
if (!identical(argsparse$wd, character(0))){
  WorkingDirectory <- toString(argsparse$wd)
}
if (!identical(argsparse$input, character(0))){
  INPUT <- toString(argsparse$input)
}
if (!identical(argsparse$samplename, character(0))){
  SAMPLENAME <- toString(argsparse$samplename)
}
if (!identical(argsparse$genome, character(0))){
  GENOME <- toString(argsparse$genome)
}

print("Arguments")
print(paste0("working directory: ", WorkingDirectory))
print(paste0("Input file: ", INPUT))
print(paste0("Sample name: ", SAMPLENAME))
print(paste0("Genome: ", GENOME))

#####################################################################
setwd(WorkingDirectory)

if (GENOME == "mm10") {
  GENOME <- "BSgenome.Mmusculus.UCSC.mm10"
  annoDb <- "org.Mm.eg.db"
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
}

if (GENOME == "hg38") {
  GENOME <- "BSgenome.Hsapiens.UCSC.hg38"
  annoDb <- "org.Hs.eg.db"
  TxDb <-TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
}

library("CAGEr")

ce <- CAGEr::CAGEexp(genomeName = GENOME, inputFiles = INPUT, inputFilesType = "bamPairedEnd", sampleLabels = SAMPLENAME)
CAGEr::getCTSS(ce,correctSystematicG=FALSE,removeFirstG=FALSE) #,useMulticore = T, nrCores = nrow(df))

gr <- CAGEr::CTSStagCountGR(ce, sample = ce@sampleMap$primary)
chosen <- sample(1:length(gr), size = 30000, replace = T, prob = gr$score %>%  as.numeric( ) )
gr <- gr[chosen,]
# here replace = T is a desired property. This command would give the effect of subsampling the reads
#instead of the regions.
gr <- GenomicRanges::sort(gr)
anno <- ChIPseeker::annotatePeak(gr, tssRegion=c(-500, 0), TxDb=TxDb, annoDb=annoDb)
# Browse[2]> t <- anno@anno %>% as.data.frame()
# Browse[2]> t$chrpos <- paste0(t$seqnames, ":", t$pos)
# Browse[2]> t$chrpos %>% duplicated() %>% which() %>% length()
# [1] 10455 # this means that the anno pipeline respects the fact that one TSS can have multiple reads

png(filename = paste0("1_", SAMPLENAME, ".TSS_distrib_pie.png"), pointsize=10, width=3200, height=2000, res=600)
ChIPseeker::plotAnnoPie(anno)
SHORT_SAMPLENAME = gsub("_BRep1_R1.fastq.gz","",gsub("Trimmed_TD_","",SAMPLENAME))
textbox(c(0,0.2), 1, SHORT_SAMPLENAME, box=FALSE, font = 0.5)
dev.off()

# save data
save.image(paste0("1_", SAMPLENAME, ".TSS_distrib_pie.RData"))
