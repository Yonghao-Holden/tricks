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
parser.add_argument("--defect", help = "running --defect output or not", type = str, default = "no")
parser.add_argument("--lim", help = "xlim and ylim", type = int, default = 400)
args = parser.parse_args()

input_path = args.input
output_file = args.output
defect = args.defect
lim = args.lim
###############################################################
# Part 1. Data Read in
###############################################################
if defect == "yes":
	input_files = glob.glob(input_path + "/" + "*lc_extrapdefects.txt")
elif defect == "no":
	input_files = glob.glob(input_path + "/" + "*lc_extrap.txt")

input_panda_list = []
for input_file in input_files:
	sample_name = os.path.splitext(os.path.basename(input_file))[0]
	input_panda = pd.read_table(input_file)
	input_panda = input_panda[['TOTAL_READS', 'EXPECTED_DISTINCT']]
	input_panda["sample"] = sample_name
	input_panda['TOTAL_READS'] = input_panda['TOTAL_READS']/1000000
	input_panda['EXPECTED_DISTINCT'] = input_panda['EXPECTED_DISTINCT']/1000000
	input_panda_list.append(input_panda)

input_panda = pd.concat(input_panda_list, ignore_index=True)

print(input_panda.head(20))
print(input_panda.shape)

###############################################################
# Part 2. Data Processing
###############################################################
fig, ax = plt.subplots(figsize=(8,8))
sns.lineplot(x = "TOTAL_READS", y = "EXPECTED_DISTINCT", hue = "sample", data = input_panda)
plt.legend('')
max_read = input_panda["TOTAL_READS"].max()
if max_read <= lim:
	plt.plot([0, max_read], [0, max_read], color = 'black', linewidth = 1, alpha=0.6, ls='--')
else:
	plt.plot([0, lim], [0, lim], color = 'black', linewidth = 1, alpha=0.6, ls='--')
plt.xlabel("Total Reads (Million)")
plt.ylabel("Distinct Reads (Million)")
plt.xlim(0,lim)
plt.ylim(0,lim)
#plt.legend(loc=3, ncol=1) #, fontsize = 6, bbox_to_anchor=(0, -1)
plt.tight_layout()
plt.savefig(output_file, dpi = 1200)
plt.close()





