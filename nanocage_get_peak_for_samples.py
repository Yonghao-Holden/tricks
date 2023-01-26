#!/usr/bin/env python3
"""
This script is to parse the output from nanocage_post_process_CAGEpeak_v2.R to obtain CTSS beak bed file for individual sample.

Usage: python3 /bar/yliang/tricks/nanocage_get_peak_for_samples.py --input ${STUDY}.capf0.15.bed --tpm 0.3 --capfilter 0.15
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
parser.add_argument("--input", help = "input file name", type = str, default = ".")
parser.add_argument("--capfilter", help = "cutoff of percentage of unannotated G reads for each peak to determine peak calling", type = float, default = "0.15")
parser.add_argument("--tpm", help = "tpm cutoff to determine peak for each sample", type = float, default = "0.3")

args = parser.parse_args()

tpm = args.tpm
input_file = args.input
capfilter = args.capfilter

###############################################################
# Part 1. Data Read in for single cell
###############################################################
input_panda = pd.read_csv(input_file, sep="\t")

for sample_index in range(7, len(input_panda.columns)):
	sample_name = input_panda.columns[sample_index]
	output_index = list(range(0,7)) + [sample_index]
	sample_panda = input_panda.iloc[:,output_index]
	sample_panda.columns = ["chr","start","end","peakname","tpm","strand","maxUGPercentage","sample_tpm"]
	sample_panda = sample_panda[["chr","start","end","peakname","sample_tpm","strand","maxUGPercentage","tpm"]]
	filtered_sample_panda = sample_panda[sample_panda["sample_tpm"] > tpm]
	filtered_sample_panda = filtered_sample_panda.astype({'start':'int', 'end':'int'})
	filtered_sample_panda.to_csv(sample_name+"_capf"+str(capfilter)+"_peak"+str(tpm)+"_bed", header=False, sep="\t", index=False)


