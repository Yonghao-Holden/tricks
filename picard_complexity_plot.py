#!/usr/bin/env python3
"""
Input folder with preseq lc_extrap result file, output plots.

Usage: python3 preseq_lc_extrap_plot.py --input --output
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

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input path", type = str)
parser.add_argument("--output", help = "output file name", type = str)
parser.add_argument("--defect", help = "running --defect output or not", type = str)
parser.add_argument("--lim", help = "xlim and ylim", type = int, default = 1000)
args = parser.parse_args()

input_path = args.input
output_file = args.output
defect = args.defect
lim = args.lim
###############################################################
# Part 1. Data Read in
###############################################################
if defect == "yes":
	input_files = glob.glob(input_path + "/" + "*.lc_extrapdefects.txt")
if defect == "no":
	input_files = glob.glob(input_path + "/" + "*.lc_extrap.txt")

input_data = {}
index = 0

for file in input_files:
	basename = os.path.splitext(os.path.basename(file))[0]
	with open(file, 'r') as file:
		file.readline()
		for line in file:
			index += 1
			entry = line.strip('\n').split('\t')
			input_data[index] = [basename, float(float(entry[0])/1000000), float(float(entry[1])/1000000)]

input_panda = pd.DataFrame.from_dict(input_data)
input_panda_T = input_panda.transpose()
input_panda_T.columns = ['sample','read','mean']
input_panda_T['read'] = input_panda_T['read'].astype(float)
input_panda_T['mean'] = input_panda_T['mean'].astype(float)
print(input_panda_T.head(20))
print(input_panda_T.shape)

###############################################################
# Part 2. Data Processing
###############################################################
## plot 3: example virtual 4C with percentage
fig, ax = plt.subplots(figsize=(8,8))
sns.lineplot(x = "read", y = "mean", hue = "sample", data = input_panda_T)
max_read = input_panda_T["read"].max()
plt.plot([0, max_read], [0, max_read], color = 'black', linewidth = 1, alpha=0.6, ls='--')
plt.xlabel("Total Reads (Million)")
plt.ylabel("Distinct Reads (Million)")
plt.xlim(0,lim)
plt.ylim(0,lim)
plt.legend(loc=3, ncol=1) #, fontsize = 6, bbox_to_anchor=(0, -1)
plt.tight_layout()
plt.savefig(output_file, dpi = 1200)
plt.close()





