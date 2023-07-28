#!/usr/bin/env python3
"""
Sequencing depth and percentage of uniquelly mapped reads

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
parser.add_argument("--output", help = "output file name", type = str, default="wgbs_SeqDepth_UniqueAlign.pdf")
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
			if "Sequence pairs analysed in total" in line:
				total_read = int(line.strip("\n").split("\t")[1])
				input_dict["number"].append(total_read/1000000)
				input_dict["data_type"].append("read_depth")
			if "Number of paired-end alignments with a unique best hit" in line:
				unique_map_read = int(line.strip("\n").split("\t")[1])
				input_dict["number"].append(unique_map_read/total_read*100)
				input_dict["data_type"].append("unique_alignment_rate")
				break

input_panda = pd.DataFrame(input_dict)
print(input_panda.head(20))
print(input_panda.shape)

###############################################################
# Part 2. Data Processing
###############################################################
input_panda = input_panda.sort_values(by='sample')
data1 = input_panda[input_panda["data_type"] == "read_depth"]
data2 = input_panda[input_panda["data_type"] == "unique_alignment_rate"]

plt.figure(figsize=(10, 6))
ax = sns.lineplot(x='sample', y='number', data=data1, color='blue', ci=None)
ax2 = ax.twinx()
sns.lineplot(x='sample', y='number', data=data2, ax=ax2, color='orange', ci=None)

ax.set_xlabel('Sample')
ax.set_ylabel('sequencing_depth (Million)', color='blue')
ax2.set_ylabel('unique_alignment_rate (%)', color='orange')

plt.title('Sequencing depth and unique alignment rate')
plt.xticks([])
plt.tight_layout()
plt.savefig(output_file, dpi = 1200)
plt.close()

input_panda.to_csv("raw_data.txt", sep="\t", index=False)



