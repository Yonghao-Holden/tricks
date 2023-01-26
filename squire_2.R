#!/usr/bin/env Rscript
library('Xmisc')
library("data.table")
library("preprocessCore")
library("EnhancedVolcano")

# This was created by Step 4 and has the data needed about each candidate.
# load('Step4.RData')

# Default parameters
WorkingDirectory <- "/scratch/yliang/HNSCC/analysis/rna-seq/1_squire/kera_vs_JHU006_NEG_call"
DATA <- "all"
OUTPUT <- "AvsB_volcano_TE.pdf"

# User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-wd',type='character',help='working directory')
parser$add_argument('-data',type='character',help='all or gene or TE')
parser$add_argument('-output',type='character',help='output file name')
argsparse <- parser$get_args()

# Argument parse to change defaults
if (!identical(argsparse$wd, character(0))){
  WorkingDirectory <- toString(argsparse$wd)
}
if (!identical(argsparse$data, character(0))){
  DATA <- toString(argsparse$data)
}
if (!identical(argsparse$output, character(0))){
  OUTPUT <- toString(argsparse$output)
}

print("Arguments")
print(paste0("working directory: ", WorkingDirectory))
print(paste0("data: ", DATA))
print(paste0("output file name: ", OUTPUT))

setwd(WorkingDirectory)

# Load 
if(data == "all"){
    input_file = "./DESeq2_all.txt"
} else if (data == "gene"){
    input_file = "./DESeq2_RefSeq_only.txt"
} else if (data == "TE"){
    input_file = "./DESeq2_TE_only.txt"
} else {
    input_file = "please enter correct data type"
}

print(input_file)
input_data = data.frame(fread(input_file))

pdf(OUTPUT)
EnhancedVolcano(input_data, lab = 'V1', x = 'log2FoldChange', y = 'pvalue')
dev.off()

# save data so you can modify the plot according to this website
# https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
save.image('volcano_plot.RData')



