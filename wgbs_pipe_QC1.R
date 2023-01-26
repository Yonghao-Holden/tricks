#!/usr/bin/env Rscript
# this script does easy analysis on wgbs data [autosome all CpG](coverage distribution, methylation distribution, pearson correlation).
library('Xmisc')
library("data.table")
library("preprocessCore")

# This was created by Step 4 and has the data needed about each candidate.
# load('Step4.RData')

# Default parameters
WorkingDirectory <- "/scratch/yliang/HNSCC/data/wgbs_BRep1/aligned"
INPUTSuffix <- "_R1_bismark_bt2_pe.deduplicated.bismark.cov.commonCpG.bed"
OUTPUT <- "wgbs_BRep1"

# User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('--wd',type='character',help='working directory')
parser$add_argument('--input',type='character',help='suffix of input files')
parser$add_argument('--output',type='character',help='study name as prefix of output file names')
argsparse <- parser$get_args()

# Argument parse to change defaults
if (!identical(argsparse$wd, character(0))){
  WorkingDirectory <- toString(argsparse$wd)
}
if (!identical(argsparse$input, character(0))){
  INPUTSuffix <- toString(argsparse$input)
}
if (!identical(argsparse$output, character(0))){
  OUTPUT <- toString(argsparse$output)
}

print("Arguments")
print(paste0("working directory: ", WorkingDirectory))
print(paste0("input file suffix: ", INPUTSuffix))
print(paste0("study name: ", OUTPUT))

setwd(WorkingDirectory)

input_files=list.files(path="./", pattern= paste0("*", INPUTSuffix), full.names=T)

input_data = list()
for (input_file in input_files){
  print(input_file)
  input_data[[input_file]] = data.frame(fread(input_file))
  }

input_matrix = Reduce(function(...)merge(...,by=c("V1","V2","V3")),input_data)

NumOfFiles = as.numeric(length(input_files))

colnames(input_matrix)[1:3] = c("chr","start","stop")
for(i in 1:NumOfFiles){
	basename = sub(paste0(".//Trimmed_WGBS_(.*)",INPUTSuffix), "\\1", input_files)[i]
	colnames(input_matrix)[(3*i+1):(3*i+3)] = c(basename, paste0(basename,"_C"), paste0(basename, "_T"))
}

methylation_matrix = input_matrix[,seq(4,length(input_matrix),3)]
library("reshape2")
melted_methylation_matrix = reshape2::melt(methylation_matrix)

#png(filename = paste0("4_", CELLLINE, "_cutoff.", CUTOFF, ".DI_delta_distribution.DMSO_consistent.zoomin_positive.png"), pointsize=10, width=3200, height=2000, res=600)
ggplot(melted_methylation_matrix, aes(x=value, color=variable)) + geom_density() + ggtitle(paste0("Methylation distribution in ", OUTPUT))
dev.off()

correlated_methylation_matrix = round(cor(methylation_matrix),2)
melted_correlated_methylation_matrix = reshape2::melt(correlated_methylation_matrix)
ggplot(data = melted_correlated_methylation_matrix, aes(x=Var1, y=Var2, fill=value)) + geom_tile()


