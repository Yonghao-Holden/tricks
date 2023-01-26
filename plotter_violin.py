#!/usr/bin/env python3
"""
Input eigenvalues in bedgraph format, output a heatmap matrix based on pearson correlation.

Usage: python3 correlation_matrix_ploting.py --input <dir> --output <output_file_basename> --corr <‘pearson’, ‘kendall’, ‘spearman’> --suffix <suffix of input files> --removenan <remove nan?>
python3 6_correlation_matrix_ploting.py --input $PWD --output B36
"""

import argparse
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import glob
import math

def mybool(s):
	return s != 'False'

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "path to input files", type = str)
parser.add_argument("--output", help = "output file name", type = str)
parser.add_argument("--col", help = "specify the column you want to plot", type = int)
parser.add_argument("--prefix", help = "prefix of input files", type = str, default = "")
parser.add_argument("--suffix", help = "suffix of input files", type = str, default = "bedgraph")
parser.add_argument("--ymax", help = "ymax", type = str, default = "false")
parser.add_argument("--ymin", help = "ymin", type = str, default = "false")

args = parser.parse_args()

input_path = args.input
output_file = args.output
col = args.col
prefix = args.prefix
suffix = args.suffix
ymax = args.ymax
ymin = args.ymin

###############################################################
# Part 1. Data Read in
###############################################################
input_files = glob.glob(input_path + "/" + prefix + "*" + suffix)

input_data = {"file":[], "data":[]}

#for file in input_files:
#	basename = os.path.splitext(os.path.basename(file))[0].split('.')[0]
#	input_data[basename] = []
#	with open(file, 'r') as file:
#		for line in file:
#			entry = line.strip('\n').split('\t')
#			input_data[basename].append(float(entry[col]))

for file in input_files:
	basename = os.path.basename(file).split(suffix)[0]
	print("%s\t%s" %(file, basename))
	with open(file, 'r') as file:
		for line in file:
			entry = line.strip('\n').split('\t')
			input_data["file"].append(basename)
			input_data["data"].append(float(entry[col-1]))

input_panda = pd.DataFrame.from_dict(input_data)

print(input_panda.head(20))
print(input_panda.shape)
###############################################################
# Part 2. Data Processing
###############################################################
plot = sns.violinplot(x = "file", y = "data", data = input_panda, cut = 0)
if ymax != "false":
	ymax = float(ymax)
	plt.ylim(top = ymax)
if ymin != "false":
	ymin = float(ymin)
	plt.ylim(bottom = ymin)
#plt.gcf().subplots_adjust(left=0.02, right=0.83, bottom = 0.17, top = 0.98)
plt.savefig(output_file, dpi=1600)
plt.close()

