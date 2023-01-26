#!/usr/bin/env python3
"""
Input fasta file with protein sequences, output potential digested peptide sequences with molecular weight information. (allowing 1, 2 missed cleavage)

Usage: python3 in_silico_digestion.py --input input.fasta
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
parser.add_argument("--input", help = "input fasta file", type = str)
parser.add_argument("--enzyme", help = "digestion enzyme: trypsin, LysC, ArgC, GluC", type = str)

args = parser.parse_args()

input_file = args.input
enzyme = args.enzyme
output_file = os.path.splitext(os.path.basename(input_file))[0] +"."+enzyme+"_digested.fa"
output_size_plot = os.path.splitext(os.path.basename(input_file))[0] +"."+enzyme+"_digested.size_dist.pdf"

def peptide_molecular_weight_calculator(seq):
	weight = 0
	for aa in seq:
		if aa=="Ala" or aa=="A":
			weight += 89.1
		if aa=="Arg" or aa=="R":
			weight += 174.2
		if aa=="Asn" or aa=="N":
			weight += 132.1
		if aa=="Asp" or aa=="D":
			weight += 133.1
		if aa=="Cys" or aa=="C":
			weight += 121.2
		if aa=="Glu" or aa=="E":
			weight += 147.1
		if aa=="Gln" or aa=="Q":
			weight += 146.2
		if aa=="Gly" or aa=="G":
			weight += 75.1
		if aa=="His" or aa=="H":
			weight += 155.2
		if aa=="Ile" or aa=="I":
			weight += 131.2
		if aa=="Leu" or aa=="L":
			weight += 131.2
		if aa=="Lys" or aa=="K":
			weight += 146.2
		if aa=="Met" or aa=="M":
			weight += 149.2
		if aa=="Phe" or aa=="F":
			weight += 165.2
		if aa=="Pro" or aa=="P":
			weight += 115.1
		if aa=="Ser" or aa=="S":
			weight += 105.1
		if aa=="Thr" or aa=="T":
			weight += 119.1
		if aa=="Trp" or aa=="W":
			weight += 204.2
		if aa=="Tyr" or aa=="Y":
			weight += 181.2
		if aa=="Val" or aa=="V":
			weight += 117.1	
	return weight

###############################################################
# Part 1. Data Read in
###############################################################
input_data = {} # input_data[index] -> [header, aa sequence]
index = 0

with open(input_file, 'r') as file:
	for line in file:
		entry = line.strip('\n')
		if entry[0] == ">":
			input_data[index] = [entry]
		else:
			input_data[index].append(entry)
			index += 1

###############################################################
# Part 2. Data Processing
# In silico digestion
# Trypsin cleaves protein peptide bonds specifically on the carboxy side of lysine (Lys, K) and arginine (Arg, R) residues
###############################################################
output_data = {}
total_peptide_count = 1

size_data = {}

for index in input_data:
	protein_name = input_data[index][0]
	protein_sequence = input_data[index][1]
	cutsite_index = [0]
	peptide_count = 1
	if enzyme == "trypsin":
		for i in range(len(protein_sequence)): 
			if protein_sequence[i] == "K" or protein_sequence[i] == "R" or protein_sequence[i] == "k" or protein_sequence[i] == "r": # try to locate K and R in the sequence
				if protein_sequence[i+1] != "P" and protein_sequence[i+1] != "p": # K/R on the left side of a P is mostly resistance to trypsin digestion
					cutsite_index.append(i)
	elif enzyme == "GluC":
		for i in range(len(protein_sequence)): 
			if protein_sequence[i] == "D" or protein_sequence[i] == "E" or protein_sequence[i] == "d" or protein_sequence[i] == "e": # try to locate D and E in the sequence
				cutsite_index.append(i)
	elif enzyme == "LysC":
		for i in range(len(protein_sequence)): 
			if protein_sequence[i] == "K" or protein_sequence[i] == "k": # try to locate K in the sequence
				cutsite_index.append(i)
	elif enzyme == "ArgC": # https://www.promega.com/products/mass-spectrometry/proteases-and-surfactants/arg_c_-sequencing-grade/?catNum=V1881
		for i in range(len(protein_sequence)): 
			if protein_sequence[i] == "K" or protein_sequence[i] == "R" or protein_sequence[i] == "k" or protein_sequence[i] == "r": # try to locate K and R in the sequence
				cutsite_index.append(i)	
	cutsite_index.append(len(protein_sequence)-1)
	# 0 missed cleavage
	for i in range(len(cutsite_index)-1):
		if i == 0:
			peptide_sequence = protein_sequence[cutsite_index[i]:cutsite_index[i+1]+1]
		else:
			peptide_sequence = protein_sequence[cutsite_index[i]+1:cutsite_index[i+1]+1]
		weight = peptide_molecular_weight_calculator(peptide_sequence)
		size_data[total_peptide_count] = [protein_name, weight]
		output_data[total_peptide_count] = [protein_name+"_peptide"+str(peptide_count)+"_"+str(cutsite_index[i]+1)+"-"+str(cutsite_index[i+1]+1)+"_"+str(len(peptide_sequence))+"_0MissedCleavage_"+str(weight)+"Da", peptide_sequence]
		peptide_count += 1
		total_peptide_count += 1
	# 1 missed cleavage
	for i in range(len(cutsite_index)-2):
		if i == 0:
			peptide_sequence = protein_sequence[cutsite_index[i]:cutsite_index[i+2]+1]
		else:
			peptide_sequence = protein_sequence[cutsite_index[i]+1:cutsite_index[i+2]+1]
		weight = peptide_molecular_weight_calculator(peptide_sequence)
		size_data[total_peptide_count] = [protein_name, weight]
		output_data[total_peptide_count] = [protein_name+"_peptide"+str(peptide_count)+"_"+str(cutsite_index[i]+1)+"-"+str(cutsite_index[i+2]+1)+"_"+str(len(peptide_sequence))+"_1MissedCleavage_"+str(weight)+"Da", peptide_sequence]
		peptide_count += 1
		total_peptide_count += 1
	# 2 missed cleavage
	for i in range(len(cutsite_index)-3):
		if i == 0:
			peptide_sequence = protein_sequence[cutsite_index[i]:cutsite_index[i+3]+1]
		else:
			peptide_sequence = protein_sequence[cutsite_index[i]+1:cutsite_index[i+3]+1]
		weight = peptide_molecular_weight_calculator(peptide_sequence)
		size_data[total_peptide_count] = [protein_name, weight]
		output_data[total_peptide_count] = [protein_name+"_peptide"+str(peptide_count)+"_"+str(cutsite_index[i]+1)+"-"+str(cutsite_index[i+3]+1)+"_"+str(len(peptide_sequence))+"_2MissedCleavage_"+str(weight)+"Da", peptide_sequence]
		peptide_count += 1
		total_peptide_count += 1

###############################################################
# Part 3. Data Output
###############################################################
with open(output_file, 'w') as file:
	for total_peptide_count in output_data:
		file.write("%s\n%s\n" %(output_data[total_peptide_count][0], output_data[total_peptide_count][1]))

size_panda = pd.DataFrame.from_dict(size_data, orient="index",columns=["sample", "size"])
sns.violinplot(x="sample", y="size", data=size_panda)
plt.axhline(y=500)
plt.axhline(y=2500)
plt.savefig(output_size_plot, dpi=1600)
plt.close()








