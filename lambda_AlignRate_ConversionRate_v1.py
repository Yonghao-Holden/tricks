#!/usr/bin/env python3
"""
Input folder with preseq lc_extrap result file, output plots.

Usage: python3 lambda_AlignRate_ConversionRate_v1.py --input . --output
"""

import argparse
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input path", type = str, default=".")
parser.add_argument("--suffix", help = "suffix of input files", type = str, default = "_bismark_bt2_PE_report.txt")
parser.add_argument("--output", help = "output file name", type = str, default="wgbs_lambda_AlignRate_ConversionRate.pdf")
args = parser.parse_args()

input_path = args.input
output_file = args.output
suffix = args.suffix

###############################################################
# Part 1. Data Read in
###############################################################
input_files = glob.glob(input_path + "/*" + suffix)

input_dict = {}
input_dict["sample"] = []
input_dict["number"] = []
input_dict["data_type"] = []
for input_file in input_files:
	sample_name = os.path.basename(input_file).rstrip(suffix)
	input_dict["sample"] += [sample_name]*2
	with open(input_file, "r") as file:
		for line in file:
			if "Mapping efficiency" in line:
				input_dict["number"].append(float(line.strip("\n").split("\t")[1].split("%")[0]))
				input_dict["data_type"].append("alignment_rate")
			if "C methylated in CpG context" in line:
				input_dict["number"].append(100-float(line.strip("\n").split("\t")[1].split("%")[0]))
				input_dict["data_type"].append("conversion_rate")
				break

input_panda = pd.DataFrame(input_dict)
print(input_panda.head(20))
print(input_panda.shape)

###############################################################
# Part 2. Data Processing
###############################################################
input_panda = input_panda.sort_values(by='sample')
data1 = input_panda[input_panda["data_type"] == "alignment_rate"]
data2 = input_panda[input_panda["data_type"] == "conversion_rate"]

plt.figure(figsize=(10, 6))
ax = sns.lineplot(x='sample', y='number', data=data1, color='blue', ci=None)
ax2 = ax.twinx()
sns.lineplot(x='sample', y='number', data=data2, ax=ax2, color='orange', ci=None)

ax.set_xlabel('Sample')
ax.set_ylabel('alignment_rate (%)', color='blue')
ax2.set_ylabel('conversion_rate (%)', color='orange')

plt.title('Alignment rate and conversion rate of lambda genome')
plt.xticks([])
plt.tight_layout()
plt.savefig(output_file, dpi = 1200)
plt.close()



