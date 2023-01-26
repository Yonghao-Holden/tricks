#!/usr/bin/env python3
"""
Get output from rna_pipe_featureCounts_to_z-score.py, plot RPKM or z-score heatmap of a set of genes

Usage: python3 rna_pipe_expression_heatmap.py --input ${STUDY}_featurecounts.z-score.txt
"""

import argparse
import sys
import os
import math
import statistics
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import re

import glob
import random as rd
from scipy import stats
import scipy

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input file", type = str, default = "rna-seq_tumor1_featurecounts.z-score.txt")
parser.add_argument("--output", help = "output file basename", type = str, default = "yo")
parser.add_argument("--sampleN", help = "sample number", type = int, default = 3)
parser.add_argument("--datatype", help = "RPKM or z-score, will be shown on plot title", type = str, default = "z-score")
parser.add_argument("--vmax", help = "max value on heatmap", type = float, default = 0)
parser.add_argument("--setting", help = "with cell type information or just a list of gene", type = str, default = "cell") # or "gene"
parser.add_argument("--genelist", help = "genelist, first column is cell type, second column is gene name", type = str, default = "/bar/yliang/genomes/private/MarkerGeneList_melted.txt") # or ", delimited list of gene names"

args = parser.parse_args()

input_file = args.input
sampleN = args.sampleN
datatype = args.datatype
genelist = args.genelist
vmax = args.vmax
setting = args.setting
output = args.output

###############################################################
# Part 0. Function
###############################################################
def get_mean(input_list):
	return sum(input_list)/len(input_list)

###############################################################
# Part 1. Data Read in
###############################################################
input_data = pd.read_csv(input_file, sep="\t", header = 0) # , index_col= 0

if datatype == "RPKM":
	samples = [input_data.columns[int(i)] for i in range(6+sampleN,6+2*sampleN,1)]	
if datatype == "log2RPKM":
	samples = [input_data.columns[int(i)] for i in range(6+2*sampleN,6+3*sampleN,1)]		
if datatype == "z-score":
	samples = [input_data.columns[int(i)] for i in range(6+3*sampleN,6+4*sampleN,1)]

simplified_samples = [re.sub("_R1\.fastqAligned\.sortedByReadname\.out\.bam_" + datatype,'',re.sub('\.\./aligned/Trimmed_mRNA_', '', sample)) for sample in samples]

if setting == "cell":
	genelist_data = {}
	with open(genelist, 'r') as file:
		header = file.readline()
		for line in file:
			entry = line.strip("\n").split("\t")
			if entry == [""]:
				continue
			if entry[1] in genelist_data:
				genelist_data[entry[1]] = genelist_data[entry[1]] + "_" + entry[0]
			if entry[1] not in genelist_data:
				genelist_data[entry[1]] = entry[0]
	flag = 1
	for gene in genelist_data.keys():
		if flag == 1:
			plot_data = input_data[input_data["Geneid"] == gene][["Geneid"]+samples]
			flag = 0
		else:
			plot_data = plot_data.append(input_data[input_data["Geneid"] == gene][["Geneid"]+samples])
	
	plot_data.columns = [re.sub("_R1\.fastqAligned\.sortedByReadname\.out\.bam_" + datatype,'',re.sub('\.\./aligned/Trimmed_mRNA_', '', sample)) for sample in plot_data.columns]
	
	plot_data["CellType"] = plot_data["Geneid"].apply(lambda x: genelist_data[x])
	
	plot_data = plot_data.set_index("Geneid")
	
	CellType = plot_data.pop("CellType")
	pallete = sns.cubehelix_palette(len(CellType.unique()),light=0.9, dark=0.1, reverse=True,start=3, rot=-4)
	lut = dict(zip(CellType.unique(), pallete))
	row_colors = CellType.map(lut)
	
	#plt.figure(figsize=(10,20))
	#cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
	if vmax == 0:
		g = sns.clustermap(plot_data, row_colors=row_colors, row_cluster=False, col_cluster=False, cmap = "vlag") #, col_colors=row_colors
	else:
		g = sns.clustermap(plot_data, row_colors=row_colors, row_cluster=False, col_cluster=False, cmap = "vlag", vmax = vmax) #, col_colors=row_colors
	
	for label in CellType.unique():
	    g.ax_col_dendrogram.bar(0, 0, color=lut[label],label=label, linewidth=0)
	g.ax_col_dendrogram.legend(loc="center", ncol=5)
	g.cax.set_position([.1, .2, .03, .45])
	#g.fig.set_figheight(20)
	#g.fig.set_figwidth(10)
	#plt.title(os.path.splitext(os.path.basename(input_file))[0] + " " + datatype)
	#plt.tight_layout()
	plt.savefig(os.path.splitext(os.path.basename(input_file))[0] + "_heatmap_" + datatype + ".pdf", dpi = 1600)
	plt.close()

if setting == "gene":
	gene_names = genelist.split(",")
	print(gene_names)
	flag = 1
	for gene_name in gene_names:
		if flag == 1:
			plot_data = input_data[input_data["Geneid"] == gene_name][["Geneid"]+samples]
			flag = 0
		else:
			plot_data = plot_data.append(input_data[input_data["Geneid"] == gene_name][["Geneid"]+samples])
	plot_data.columns = [re.sub("_R1\.fastqAligned\.sortedByReadname\.out\.bam_" + datatype,'',re.sub('\.\./aligned/Trimmed_mRNA_', '', sample)) for sample in plot_data.columns]
	plot_data = plot_data.set_index("Geneid")
	print(plot_data)
	if vmax == 0:
		g = sns.clustermap(plot_data, row_cluster=False, col_cluster=False, cmap = "vlag", cbar_kws={'label': datatype}) #, col_colors=row_colors
	else:
		g = sns.clustermap(plot_data, row_cluster=False, col_cluster=False, cmap = "vlag", cbar_kws={'label': datatype}, vmax = vmax) #, col_colors=row_colors
	if output == "yo":
		plt.savefig(os.path.splitext(os.path.basename(input_file))[0] + "_heatmap_" + datatype + "_" + '-'.join(gene_names) + ".pdf", dpi = 1600)
	else:
		plt.savefig(os.path.splitext(os.path.basename(input_file))[0] + "_heatmap_" + datatype + "_" + output + ".pdf", dpi = 1600)
	plt.close()












































