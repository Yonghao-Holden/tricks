#!/usr/bin/env python3
"""
Convert a table into a heatmap

Usage: python3 heatmap.py
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
import random as rd
import math
from scipy import stats
import scipy
from matplotlib.colors import LinearSegmentedColormap

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input file", type = str)
parser.add_argument("--output", help = "output file", type = str)
parser.add_argument("--color", help = "vmax color", type = str, default = "default")
parser.add_argument("--center", help = "center color", type = float, default = 999)
parser.add_argument("--vmin", help = "vmin for heatmap", type = float, default = 999)
parser.add_argument("--vmax", help = "vmax for heatmap", type = float, default = 999)
parser.add_argument("--label", help = "color bar label", type = str)
parser.add_argument("--title", help = "title", type = str)

args = parser.parse_args()

input_file = args.input
output_file = args.output
color = args.color
center = args.center
vmin = args.vmin
vmax = args.vmax
label = args.label
title =args.title

###############################################################
# Part 1. Data Read in
###############################################################
data = {}

with open(input_file, 'r') as file:
	header = file.readline().strip('\n').split('\t')[1:]
	for line in file:
		entry = line.strip('\n').split()
		data[entry[0]] = [float(x) for x in entry[1:]]

data_panda = pd.DataFrame.from_dict(data).transpose()
data_panda.columns = header
print(data_panda)

#cmap = sns.dark_palette("#4C72B0", as_cmap=True)
plt.figure(figsize=(5,6))
if color != "default": # custom color
	if center != 999: # have center value
		cmap = LinearSegmentedColormap.from_list(name = "test", colors=['black', 'white', color])
		if vmin != 999 or vmax != 999: # have custom vmin and vmax
			sns.heatmap(data_panda, cmap = cmap, center = center, vmin = vmin, vmax = vmax, cbar_kws={'label': label})
		else:
			sns.heatmap(data_panda, cmap = cmap, center = center, cbar_kws={'label': label})
	elif center == 999: # don't have center value
		cmap = LinearSegmentedColormap.from_list(name = "test", colors=['black', color])
		if vmin != 999 or vmax != 999: # have custom vmin and vmax		
			sns.heatmap(data_panda, cmap = cmap, vmin = vmin, vmax = vmax, cbar_kws={'label': label})
		else:
			sns.heatmap(data_panda, cmap = cmap, cbar_kws={'label': label})
elif color == "default":
	if center != 999:
		sns.heatmap(data_panda, center = center, cbar_kws={'label': label})
	elif center == 999:
		sns.heatmap(data_panda, cbar_kws={'label': label})
if title != "yo":
	plt.title(title)
plt.tight_layout()
plt.savefig(output_file, dpi = 1200)
plt.close()


