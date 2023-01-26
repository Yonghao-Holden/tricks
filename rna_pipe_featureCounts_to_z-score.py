#!/usr/bin/env python3
"""
Get output from featureCounts, calculate RPKM and z-score.

Usage: python3 featureCounts_to_z-score.py --input ${STUDY}_featurecounts.txt
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
from sklearn.decomposition import PCA
import math


import glob
import random as rd
from scipy import stats
import scipy

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input file", type = str)
parser.add_argument("--size", help = "file with library size information", type = str)
parser.add_argument("--PC", help = "dimension to be plotted", type = str, default = "1,2")

args = parser.parse_args()

input_file = args.input
size_file = args.size
dimension_num = [int(x) for x in args.PC.split(',')]
dimension = ["PC" + x for x in args.PC.split(',')]

###############################################################
# Part 0. Function
###############################################################
def get_mean(input_list):
	return sum(input_list)/len(input_list)

###############################################################
# Part 1. Data Read in
###############################################################
input_data = pd.read_csv(input_file, sep="\t", header = 1, skiprows = 0) #, index_col= 0

sample_size = {}
with open(size_file, 'r') as file:
	for line in file:
		entry = line.strip('\n').split('\t')
		sample_size[entry[0]] = int(entry[1])

sample_num = len(input_data.columns) - 6

for i in range(sample_num):
	sample_name =  input_data.columns[i+6]
	#total_read_num = sum(input_data.iloc[:,i+6])/1000000 # total number of reads that overlap with exons
	total_read_num = sample_size[sample_name] # total number of uniquely mapped reads in the library (bam file)
	input_data[sample_name + "_RPKM"] = (input_data[sample_name] / total_read_num) / (input_data["Length"]/1000) # RPKM

for i in range(sample_num):
	sample_name =  input_data.columns[i+6]	
	input_data[sample_name + "_log2RPKM"] = np.log2(input_data[sample_name + "_RPKM"]+1)

for index, gene in input_data.iterrows():
	gene_data = list(input_data.iloc[index,6+2*(sample_num):6+3*(sample_num)])
	std = statistics.stdev(gene_data)
	mean = get_mean(gene_data)
	for i in range(sample_num):
		sample_name =  input_data.columns[i+6]
		if std != 0:
			input_data.loc[index, sample_name + "_z-score"] = (gene_data[i] - mean)/std
		if std == 0:
			input_data.loc[index, sample_name + "_z-score"] = 0


input_data.to_csv(os.path.splitext(os.path.basename(input_file))[0]+".z-score.txt", sep="\t", index = False)

###############################################################
# Part 2. correlation heatmap (use RPKM as input)
###############################################################
plot_data = input_data.iloc[:,6+sample_num:6+2*sample_num] ## RPKM

plot_data.columns = [re.sub("_R1\.fastqAligned\.sortedByReadname\.out\.bam_RPKM",'',re.sub('\.\./aligned/Trimmed_mRNA_', '', sample)) for sample in plot_data.columns]

sns.clustermap(plot_data.corr(), cmap = "vlag")
plt.savefig(os.path.splitext(os.path.basename(input_file))[0] + "_CorrHeatmap.pdf", dpi = 1600)
plt.close()

###############################################################
# Part 3. PCA plot (Use RPKM as input)
###############################################################
sns.set(color_codes=True)
pca = PCA()
pca.fit(plot_data.T)
pca_data = pca.transform(plot_data.T)

per_var = np.round(pca.explained_variance_ratio_*100, decimals = 1)

labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]

plt.bar(x=range(1,len(per_var)+1), height=per_var, tick_label=labels)
plt.ylabel('Percentage of Explained Variantce')
plt.xlabel('Principle Component')
plt.title('Scree Plot')
plt.tight_layout()
plt.savefig(os.path.splitext(os.path.basename(input_file))[0] + "_scree.pdf", dpi = 1600)
plt.close()

pca_df = pd.DataFrame(pca_data, index = plot_data.columns, columns = labels)
 
print(pca_df)

sns.scatterplot(pca_df.loc[:, dimension[0]], pca_df.loc[:, dimension[1]], alpha = 0.6)
plt.title(os.path.splitext(os.path.basename(input_file))[0] + " PCA plot")
plt.xlabel(dimension[0] + ' - {0}%'.format(per_var[dimension_num[0]-1]))
plt.ylabel(dimension[1] + ' - {0}%'.format(per_var[dimension_num[1]-1]))
for sample in pca_df.index:
	plt.annotate(sample, (pca_df.loc[sample,dimension[0]], pca_df.loc[sample,dimension[1]]), size = 6)
plt.legend(loc=1, ncol=1, fontsize = 6, bbox_to_anchor=(1.25, 1))
plt.tight_layout()
plt.savefig(os.path.splitext(os.path.basename(input_file))[0] + "_PCA.pdf", dpi = 1600)
plt.close()



