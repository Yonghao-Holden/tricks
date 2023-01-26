#!/usr/bin/env python3
"""
Run it under the folder as shown below, output a xlsx file with all the data summarized from inter.txt and inter_30.txt

--> run it here
	--> cell 1
		--> aligned
		--> debug
	--> cell 2
		--> aligned
		--> debug

Usage: python3 hic_qc_summary.py --input $PWD --output $STUDY
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
parser.add_argument("--input", help = "path to input files", type = str, default = ".")
parser.add_argument("--prefix", help = "prefix of input files", type = str, default = "")
parser.add_argument("--output", help = "output file basename", type = str)
parser.add_argument("--datatype", help = "single or bulk", type = str, default = "single") # options are: single, bulk, bulk_mega

args = parser.parse_args()

prefix = args.prefix
input_path = args.input
output_file = args.output
datatype = args.datatype

###############################################################
# Part 1. Data Read in for single cell
###############################################################
if(datatype == "single"):
	input_files = glob.glob(input_path + "/" + prefix + "*/aligned/inter.txt")
	input_files_30 = glob.glob(input_path + "/" + prefix + "*/aligned/inter_30.txt")
	
	input_data = {}
	input_data_30 = {}
	colnames = ["sample_name", "Sequenced Read Pairs", "Normal Paired", "Normal Paired (%)", "Chimeric Paired", "Chimeric Paired (%)", "Chimeric Ambiguous", "Chimeric Ambiguous (%)", "Unmapped", "Unmapped (%)", "Ligation Motif Present", "Ligation Motif Present (%)", "Alignable (Normal+Chimeric Paired)", "Alignable (Normal+Chimeric Paired) (%)", "Final # of interaction used for hic file generation", "Final # of interaction used for hic file generation (%/total)", "WGA Duplicates", "WGA Duplicates (%)", "Optical Duplicates", "Optical Duplicates (%)", "Library Complexity Estimate", "Intra-fragment Reads", "Intra-fragment Reads (%/total)", "Intra-fragment Reads (%/final)", "Below MAPQ Threshold", "Below MAPQ Threshold (%/total)", "Below MAPQ Threshold (%/final)", "Hi-C Contacts", "Hi-C Contacts (%/total)", "Hi-C Contacts (%/final)", "Ligation Motif Present", "Ligation Motif Present (%/total)", "Ligation Motif Present (%/final)", "3' Bias (Long Range)", "Pair Type %(L-I-O-R)", "Inter-chromosomal", "Inter-chromosomal (%/total)", "Inter-chromosomal (%/final)", "Intra-chromosomal", "Intra-chromosomal (%/total)", "Intra-chromosomal (%/final)", "Short Range (<20Kb)", "Short Range (<20Kb) (%/total)", "Short Range (<20Kb) (%/final)", "Long Range (>20Kb)", "Long Range (>20Kb) (%/total)", "Long Range (>20Kb) (%/final)", "PCR Duplicates", "Read filtering for<=8 interactions per bin"]
	
	for input_file in input_files:
		sample_name = input_file.split('/')[-3]
		if sample_name == "mega":
			continue
		input_data[sample_name] = [sample_name]
		with open(input_file, 'r') as file:
			line_counter = 0
			for a in range(1):
				b = file.readline()
			for line in file:
				line_counter+=1
				entry = line.strip('\n').split(':')[1].split()
				if line_counter in [1, 11, 22, 23, 24]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
				elif line_counter in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
					input_data[sample_name].append(entry[1].split('(')[1].split(')')[0])
				elif line_counter in [12,13,14,15,18,19,20,21]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
					input_data[sample_name].append(entry[1].split('(')[1])
					input_data[sample_name].append(entry[3].split(')')[0])
				elif line_counter in [16, 17]:
					entry = line.strip('\n').split(':')[1]
					input_data[sample_name].append(entry)
	
	for input_file in input_files_30:
		sample_name = input_file.split('/')[-3]
		if sample_name == "mega":
			continue
		input_data_30[sample_name] = [sample_name]
		with open(input_file, 'r') as file:
			line_counter = 0
			for a in range(1):
				b = file.readline()
			for line in file:
				line_counter+=1
				entry = line.strip('\n').split(':')[1].split()
				if line_counter in [1, 11, 22, 23, 24]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
				elif line_counter in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
					input_data_30[sample_name].append(entry[1].split('(')[1].split(')')[0])
				elif line_counter in [12,13,14,15,18,19,20,21]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
					input_data_30[sample_name].append(entry[1].split('(')[1])
					input_data_30[sample_name].append(entry[3].split(')')[0])
				elif line_counter in [16, 17]:
					entry = line.strip('\n').split(':')[1]
					input_data_30[sample_name].append(entry)
	
	#with open(input_file, 'r') as file:
	#	for a in range(1):
	#		b = file.readline()
	#	for line in file:
	#		entry = line.strip('\n').split(':')
	#		colnames.append(str(entry[0]))
	
	input_panda = pd.DataFrame.from_dict(input_data).T
	input_panda_30 = pd.DataFrame.from_dict(input_data_30).T
	
	input_panda.columns = colnames
	input_panda_30.columns = colnames
	
	input_panda = input_panda.sort_index()
	input_panda_30 = input_panda_30.sort_index()
	
	input_panda = input_panda.rename(columns = {'WGA Duplicates':'All Duplicates'})
	input_panda_30 = input_panda_30.rename(columns = {'WGA Duplicates':'All Duplicates'})
	input_panda = input_panda.rename(columns = {'WGA Duplicates (%)':'All Duplicates (%/total)'})
	input_panda_30 = input_panda_30.rename(columns = {'WGA Duplicates (%)':'All Duplicates (%/total)'})
	
	input_panda["Unique Reads"] = input_panda["Alignable (Normal+Chimeric Paired)"] - input_panda["All Duplicates"] - input_panda["Optical Duplicates"]
	input_panda["Unique Reads (%/total)"] = input_panda["Unique Reads"]/input_panda["Sequenced Read Pairs"]*100
	input_panda["PCR Duplicates (%/total)"] = input_panda["PCR Duplicates"]/input_panda["Sequenced Read Pairs"]*100
	input_panda["PCR Duplicates (%/alignable)"] = input_panda["PCR Duplicates"]/input_panda["Alignable (Normal+Chimeric Paired)"]*100
	input_panda["Final # of interaction used for hic file generation (%/unique)"] = input_panda["Final # of interaction used for hic file generation"]/input_panda["Unique Reads"]*100
	input_panda["Read filtering for<=8 interactions per bin (%/total)"] = input_panda["Read filtering for<=8 interactions per bin"]/input_panda["Sequenced Read Pairs"]*100
	input_panda["Read filtering for<=8 interactions per bin (%/unique)"] = input_panda["Read filtering for<=8 interactions per bin"]/input_panda["Unique Reads"]*100
	input_panda["All Duplicates (%/alignable)"] = input_panda["All Duplicates"]/input_panda["Alignable (Normal+Chimeric Paired)"]*100
	input_panda["WGA Duplicates"] = input_panda["All Duplicates"] - input_panda["PCR Duplicates"]
	input_panda["WGA Duplicates (%/total)"] = input_panda["WGA Duplicates"]/input_panda["Sequenced Read Pairs"]*100
	input_panda["WGA Duplicates (%/alignable)"] = input_panda["WGA Duplicates"]/input_panda["Alignable (Normal+Chimeric Paired)"]*100
	input_panda["Short Range (<20Kb) (%/Hi-C contacts)"] = input_panda["Short Range (<20Kb)"]/input_panda["Hi-C Contacts"]*100
	input_panda["Long Range (>20Kb) (%/Hi-C contacts)"] = input_panda["Long Range (>20Kb)"]/input_panda["Hi-C Contacts"]*100
	input_panda["Inter-chromosomal (%/Hi-C contacts)"] = input_panda["Inter-chromosomal"]/input_panda["Hi-C Contacts"]*100
	
	input_panda_30["Unique Reads"] = input_panda_30["Alignable (Normal+Chimeric Paired)"] - input_panda_30["All Duplicates"] - input_panda_30["Optical Duplicates"]
	input_panda_30["Unique Reads (%/total)"] = input_panda_30["Unique Reads"]/input_panda_30["Sequenced Read Pairs"]*100
	input_panda_30["PCR Duplicates (%/total)"] = input_panda_30["PCR Duplicates"]/input_panda_30["Sequenced Read Pairs"]*100
	input_panda_30["PCR Duplicates (%/alignable)"] = input_panda_30["PCR Duplicates"]/input_panda_30["Alignable (Normal+Chimeric Paired)"]*100
	input_panda_30["Final # of interaction used for hic file generation (%/unique)"] = input_panda_30["Final # of interaction used for hic file generation"]/input_panda_30["Unique Reads"]*100
	input_panda_30["Read filtering for<=8 interactions per bin (%/total)"] = input_panda_30["Read filtering for<=8 interactions per bin"]/input_panda_30["Sequenced Read Pairs"]*100
	input_panda_30["Read filtering for<=8 interactions per bin (%/unique)"] = input_panda_30["Read filtering for<=8 interactions per bin"]/input_panda_30["Unique Reads"]*100
	input_panda_30["All Duplicates (%/alignable)"] = input_panda_30["All Duplicates"]/input_panda_30["Alignable (Normal+Chimeric Paired)"]*100
	input_panda_30["WGA Duplicates"] = input_panda_30["All Duplicates"] - input_panda_30["PCR Duplicates"]
	input_panda_30["WGA Duplicates (%/total)"] = input_panda_30["WGA Duplicates"]/input_panda_30["Sequenced Read Pairs"]*100
	input_panda_30["WGA Duplicates (%/alignable)"] = input_panda_30["WGA Duplicates"]/input_panda_30["Alignable (Normal+Chimeric Paired)"]*100
	input_panda_30["Short Range (<20Kb) (%/Hi-C contacts)"] = input_panda_30["Short Range (<20Kb)"]/input_panda_30["Hi-C Contacts"]*100
	input_panda_30["Long Range (>20Kb) (%/Hi-C contacts)"] = input_panda_30["Long Range (>20Kb)"]/input_panda_30["Hi-C Contacts"]*100
	input_panda_30["Inter-chromosomal (%/Hi-C contacts)"] = input_panda_30["Inter-chromosomal"]/input_panda_30["Hi-C Contacts"]*100
	
	input_panda["Unique Reads"] = input_panda["Unique Reads"].astype(float)
	input_panda["Unique Reads (%/total)"] = input_panda["Unique Reads (%/total)"].astype(float).round(decimals=2)
	input_panda["PCR Duplicates (%/total)"] = input_panda["PCR Duplicates (%/total)"].astype(float).round(decimals=2)
	input_panda["PCR Duplicates (%/alignable)"] = input_panda["PCR Duplicates (%/alignable)"].astype(float).round(decimals=2)
	input_panda["Final # of interaction used for hic file generation (%/unique)"] = input_panda["Final # of interaction used for hic file generation (%/unique)"].astype(float).round(decimals=2)
	input_panda["Read filtering for<=8 interactions per bin (%/total)"] = input_panda["Read filtering for<=8 interactions per bin (%/total)"].astype(float).round(decimals=2)
	input_panda["Read filtering for<=8 interactions per bin (%/unique)"] = input_panda["Read filtering for<=8 interactions per bin (%/unique)"].astype(float).round(decimals=2)
	input_panda["All Duplicates (%/alignable)"] = input_panda["All Duplicates (%/alignable)"].astype(float).round(decimals=2)
	input_panda["WGA Duplicates"] = input_panda["WGA Duplicates"].astype(float)
	input_panda["WGA Duplicates (%/total)"] = input_panda["WGA Duplicates (%/total)"].astype(float).round(decimals=2)
	input_panda["WGA Duplicates (%/alignable)"] = input_panda["WGA Duplicates (%/alignable)"].astype(float).round(decimals=2)
	input_panda["Short Range (<20Kb) (%/Hi-C contacts)"] = input_panda["Short Range (<20Kb) (%/Hi-C contacts)"].astype(float).round(decimals=2)
	input_panda["Long Range (>20Kb) (%/Hi-C contacts)"] = input_panda["Long Range (>20Kb) (%/Hi-C contacts)"].astype(float).round(decimals=2)
	input_panda["Inter-chromosomal (%/Hi-C contacts)"] = input_panda["Inter-chromosomal (%/Hi-C contacts)"].astype(float).round(decimals=2)
	
	input_panda_30["Unique Reads"] = input_panda_30["Unique Reads"].astype(float)
	input_panda_30["Unique Reads (%/total)"] = input_panda_30["Unique Reads (%/total)"].astype(float).round(decimals=2)
	input_panda_30["PCR Duplicates (%/total)"] = input_panda_30["PCR Duplicates (%/total)"].astype(float).round(decimals=2)
	input_panda_30["PCR Duplicates (%/alignable)"] = input_panda_30["PCR Duplicates (%/alignable)"].astype(float).round(decimals=2)
	input_panda_30["Final # of interaction used for hic file generation (%/unique)"] = input_panda_30["Final # of interaction used for hic file generation (%/unique)"].astype(float).round(decimals=2)
	input_panda_30["Read filtering for<=8 interactions per bin (%/total)"] = input_panda_30["Read filtering for<=8 interactions per bin (%/total)"].astype(float).round(decimals=2)
	input_panda_30["Read filtering for<=8 interactions per bin (%/unique)"] = input_panda_30["Read filtering for<=8 interactions per bin (%/unique)"].astype(float).round(decimals=2)
	input_panda_30["All Duplicates (%/alignable)"] = input_panda_30["All Duplicates (%/alignable)"].astype(float).round(decimals=2)
	input_panda_30["WGA Duplicates"] = input_panda_30["WGA Duplicates"].astype(float)
	input_panda_30["WGA Duplicates (%/total)"] = input_panda_30["WGA Duplicates (%/total)"].astype(float).round(decimals=2)
	input_panda_30["WGA Duplicates (%/alignable)"] = input_panda_30["WGA Duplicates (%/alignable)"].astype(float).round(decimals=2)
	input_panda_30["Short Range (<20Kb) (%/Hi-C contacts)"] = input_panda_30["Short Range (<20Kb) (%/Hi-C contacts)"].astype(float).round(decimals=2)
	input_panda_30["Long Range (>20Kb) (%/Hi-C contacts)"] = input_panda_30["Long Range (>20Kb) (%/Hi-C contacts)"].astype(float).round(decimals=2)
	input_panda_30["Inter-chromosomal (%/Hi-C contacts)"] = input_panda_30["Inter-chromosomal (%/Hi-C contacts)"].astype(float).round(decimals=2)
	
	input_panda["Unique Reads (%/total)"] = input_panda["Unique Reads (%/total)"].astype(str) + "%"
	input_panda["PCR Duplicates (%/total)"] = input_panda["PCR Duplicates (%/total)"].astype(str) + "%"
	input_panda["PCR Duplicates (%/alignable)"] = input_panda["PCR Duplicates (%/alignable)"].astype(str) + "%"
	input_panda["Final # of interaction used for hic file generation (%/unique)"] = input_panda["Final # of interaction used for hic file generation (%/unique)"].astype(str) + "%"
	input_panda["Read filtering for<=8 interactions per bin (%/total)"] = input_panda["Read filtering for<=8 interactions per bin (%/total)"].astype(str) + "%"
	input_panda["Read filtering for<=8 interactions per bin (%/unique)"] = input_panda["Read filtering for<=8 interactions per bin (%/unique)"].astype(str) + "%"
	input_panda["All Duplicates (%/alignable)"] = input_panda["All Duplicates (%/alignable)"].astype(str) + "%"
	input_panda["WGA Duplicates (%/total)"] = input_panda["WGA Duplicates (%/total)"].astype(str) + "%"
	input_panda["WGA Duplicates (%/alignable)"] = input_panda["WGA Duplicates (%/alignable)"].astype(str) + "%"
	input_panda["Short Range (<20Kb) (%/Hi-C contacts)"] = input_panda["Short Range (<20Kb) (%/Hi-C contacts)"].astype(str) + "%"
	input_panda["Long Range (>20Kb) (%/Hi-C contacts)"] = input_panda["Long Range (>20Kb) (%/Hi-C contacts)"].astype(str) + "%"
	input_panda["Inter-chromosomal (%/Hi-C contacts)"] = input_panda["Inter-chromosomal (%/Hi-C contacts)"].astype(str) + "%"
	
	input_panda_30["Unique Reads (%/total)"] = input_panda_30["Unique Reads (%/total)"].astype(str) + "%"
	input_panda_30["PCR Duplicates (%/total)"] = input_panda_30["PCR Duplicates (%/total)"].astype(str) + "%"
	input_panda_30["PCR Duplicates (%/alignable)"] = input_panda_30["PCR Duplicates (%/alignable)"].astype(str) + "%"
	input_panda_30["Final # of interaction used for hic file generation (%/unique)"] = input_panda_30["Final # of interaction used for hic file generation (%/unique)"].astype(str) + "%"
	input_panda_30["Read filtering for<=8 interactions per bin (%/total)"] = input_panda_30["Read filtering for<=8 interactions per bin (%/total)"].astype(str) + "%"
	input_panda_30["Read filtering for<=8 interactions per bin (%/unique)"] = input_panda_30["Read filtering for<=8 interactions per bin (%/unique)"].astype(str) + "%"
	input_panda_30["All Duplicates (%/alignable)"] = input_panda_30["All Duplicates (%/alignable)"].astype(str) + "%"
	input_panda_30["WGA Duplicates (%/total)"] = input_panda_30["WGA Duplicates (%/total)"].astype(str) + "%"
	input_panda_30["WGA Duplicates (%/alignable)"] = input_panda_30["WGA Duplicates (%/alignable)"].astype(str) + "%"
	input_panda_30["Short Range (<20Kb) (%/Hi-C contacts)"] = input_panda_30["Short Range (<20Kb) (%/Hi-C contacts)"].astype(str) + "%"
	input_panda_30["Long Range (>20Kb) (%/Hi-C contacts)"] = input_panda_30["Long Range (>20Kb) (%/Hi-C contacts)"].astype(str) + "%"
	input_panda_30["Inter-chromosomal (%/Hi-C contacts)"] = input_panda_30["Inter-chromosomal (%/Hi-C contacts)"].astype(str) + "%"
	
	input_panda = input_panda[["sample_name", "Sequenced Read Pairs", "Normal Paired", "Normal Paired (%)", "Chimeric Paired", "Chimeric Paired (%)", "Chimeric Ambiguous", "Chimeric Ambiguous (%)", "Unmapped", "Unmapped (%)", "Ligation Motif Present", "Ligation Motif Present (%)", "Alignable (Normal+Chimeric Paired)", "Alignable (Normal+Chimeric Paired) (%)", "PCR Duplicates", "PCR Duplicates (%/total)", "PCR Duplicates (%/alignable)", "WGA Duplicates", "WGA Duplicates (%/total)", "WGA Duplicates (%/alignable)", 'All Duplicates', "All Duplicates (%/total)", "All Duplicates (%/alignable)", "Optical Duplicates", "Optical Duplicates (%)", "Unique Reads", "Unique Reads (%/total)", "Library Complexity Estimate", "Read filtering for<=8 interactions per bin", "Read filtering for<=8 interactions per bin (%/total)", "Read filtering for<=8 interactions per bin (%/unique)", "Final # of interaction used for hic file generation", "Final # of interaction used for hic file generation (%/total)", "Final # of interaction used for hic file generation (%/unique)", "Intra-fragment Reads", "Intra-fragment Reads (%/total)", "Intra-fragment Reads (%/final)", "Below MAPQ Threshold", "Below MAPQ Threshold (%/total)", "Below MAPQ Threshold (%/final)", "Hi-C Contacts", "Hi-C Contacts (%/total)", "Hi-C Contacts (%/final)", "Ligation Motif Present", "Ligation Motif Present (%/total)", "Ligation Motif Present (%/final)", "3' Bias (Long Range)", "Pair Type %(L-I-O-R)", "Inter-chromosomal", "Inter-chromosomal (%/total)", "Inter-chromosomal (%/final)", "Inter-chromosomal (%/Hi-C contacts)", "Intra-chromosomal", "Intra-chromosomal (%/total)", "Intra-chromosomal (%/final)", "Short Range (<20Kb)", "Short Range (<20Kb) (%/total)", "Short Range (<20Kb) (%/final)", "Short Range (<20Kb) (%/Hi-C contacts)", "Long Range (>20Kb)", "Long Range (>20Kb) (%/total)", "Long Range (>20Kb) (%/final)", "Long Range (>20Kb) (%/Hi-C contacts)"]]
	input_panda_30 = input_panda_30[["sample_name", "Sequenced Read Pairs", "Normal Paired", "Normal Paired (%)", "Chimeric Paired", "Chimeric Paired (%)", "Chimeric Ambiguous", "Chimeric Ambiguous (%)", "Unmapped", "Unmapped (%)", "Ligation Motif Present", "Ligation Motif Present (%)", "Alignable (Normal+Chimeric Paired)", "Alignable (Normal+Chimeric Paired) (%)", "PCR Duplicates", "PCR Duplicates (%/total)", "PCR Duplicates (%/alignable)", "WGA Duplicates", "WGA Duplicates (%/total)", "WGA Duplicates (%/alignable)", 'All Duplicates', "All Duplicates (%/total)", "All Duplicates (%/alignable)", "Optical Duplicates", "Optical Duplicates (%)", "Unique Reads", "Unique Reads (%/total)", "Library Complexity Estimate", "Read filtering for<=8 interactions per bin", "Read filtering for<=8 interactions per bin (%/total)", "Read filtering for<=8 interactions per bin (%/unique)", "Final # of interaction used for hic file generation", "Final # of interaction used for hic file generation (%/total)", "Final # of interaction used for hic file generation (%/unique)", "Intra-fragment Reads", "Intra-fragment Reads (%/total)", "Intra-fragment Reads (%/final)", "Below MAPQ Threshold", "Below MAPQ Threshold (%/total)", "Below MAPQ Threshold (%/final)", "Hi-C Contacts", "Hi-C Contacts (%/total)", "Hi-C Contacts (%/final)", "Ligation Motif Present", "Ligation Motif Present (%/total)", "Ligation Motif Present (%/final)", "3' Bias (Long Range)", "Pair Type %(L-I-O-R)", "Inter-chromosomal", "Inter-chromosomal (%/total)", "Inter-chromosomal (%/final)", "Inter-chromosomal (%/Hi-C contacts)", "Intra-chromosomal", "Intra-chromosomal (%/total)", "Intra-chromosomal (%/final)", "Short Range (<20Kb)", "Short Range (<20Kb) (%/total)", "Short Range (<20Kb) (%/final)", "Short Range (<20Kb) (%/Hi-C contacts)", "Long Range (>20Kb)", "Long Range (>20Kb) (%/total)", "Long Range (>20Kb) (%/final)", "Long Range (>20Kb) (%/Hi-C contacts)"]]
	
	with pd.ExcelWriter(output_file+'_qc.xlsx') as writer:  
	    input_panda.to_excel(writer, sheet_name='inter')
	    input_panda_30.to_excel(writer, sheet_name='inter_30')
	
	input_panda.to_csv(output_file+'_inter.tsv', sep='\t')
	input_panda_30.to_csv(output_file+'_inter_30.tsv', sep='\t')

if(datatype == "bulk"):
	input_files = glob.glob(input_path + "/" + prefix + "*/aligned/inter.txt")
	input_files_30 = glob.glob(input_path + "/" + prefix + "*/aligned/inter_30.txt")
	
	input_data = {}
	input_data_30 = {}
	colnames = ["sample_name", "Sequenced Read Pairs", "Normal Paired", "Normal Paired (%)", "Chimeric Paired", "Chimeric Paired (%)", "Chimeric Ambiguous", "Chimeric Ambiguous (%)", "Unmapped", "Unmapped (%)", "Ligation Motif Present", "Ligation Motif Present (%)", "Alignable (Normal+Chimeric Paired)", "Alignable (Normal+Chimeric Paired) (%)", "Unique Reads", "Unique Reads (%/total)",  "PCR Duplicates",  "PCR Duplicates (%/total)", "Optical Duplicates", "Optical Duplicates (%/total)", "Library Complexity Estimate", "Intra-fragment Reads", "Intra-fragment Reads (%/total)", "Intra-fragment Reads (%/unique)", "Below MAPQ Threshold", "Below MAPQ Threshold (%/total)", "Below MAPQ Threshold (%/unique)", "Hi-C Contacts", "Hi-C Contacts (%/total)", "Hi-C Contacts (%/unique)", "Ligation Motif Present", "Ligation Motif Present (%/total)", "Ligation Motif Present (%/unique)", "3' Bias (Long Range)", "Pair Type %(L-I-O-R)", "Inter-chromosomal", "Inter-chromosomal (%/total)", "Inter-chromosomal (%/unique)", "Intra-chromosomal", "Intra-chromosomal (%/total)", "Intra-chromosomal (%/unique)", "Short Range (<20Kb)", "Short Range (<20Kb) (%/total)", "Short Range (<20Kb) (%/unique)", "Long Range (>20Kb)", "Long Range (>20Kb) (%/total)", "Long Range (>20Kb) (%/unique)"]
	
	for input_file in input_files:
		sample_name = input_file.split('/')[-3]
		if sample_name == "mega":
			continue
		input_data[sample_name] = [sample_name]
		with open(input_file, 'r') as file:
			line_counter = 0
			for a in range(1):
				b = file.readline() # remove header
			for line in file:
				line_counter+=1
				entry = line.strip('\n').split(':')[1].split()
				if line_counter in [1, 11]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
				elif line_counter in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
					input_data[sample_name].append(entry[1].split('(')[1].split(')')[0])
				elif line_counter in [12,13,14,15,18,19,20,21]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
					input_data[sample_name].append(entry[1].split('(')[1])
					input_data[sample_name].append(entry[3].split(')')[0])
				elif line_counter in [16, 17]:
					entry = line.strip('\n').split(':')[1]
					input_data[sample_name].append(entry)
	
	for input_file in input_files_30:
		sample_name = input_file.split('/')[-3]
		if sample_name == "mega":
			continue
		input_data_30[sample_name] = [sample_name]
		with open(input_file, 'r') as file:
			line_counter = 0
			for a in range(1):
				b = file.readline() # remove header
			for line in file:
				line_counter+=1
				entry = line.strip('\n').split(':')[1].split()
				if line_counter in [1, 11]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
				elif line_counter in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
					input_data_30[sample_name].append(entry[1].split('(')[1].split(')')[0])
				elif line_counter in [12,13,14,15,18,19,20,21]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
					input_data_30[sample_name].append(entry[1].split('(')[1])
					input_data_30[sample_name].append(entry[3].split(')')[0])
				elif line_counter in [16, 17]:
					entry = line.strip('\n').split(':')[1]
					input_data_30[sample_name].append(entry)
	
	input_panda = pd.DataFrame.from_dict(input_data).T
	input_panda_30 = pd.DataFrame.from_dict(input_data_30).T
	
	input_panda.columns = colnames
	input_panda_30.columns = colnames
	
	input_panda = input_panda.sort_index()
	input_panda_30 = input_panda_30.sort_index()
	
	with pd.ExcelWriter(output_file+'_qc.xlsx') as writer:  
	    input_panda.to_excel(writer, sheet_name='inter')
	    input_panda_30.to_excel(writer, sheet_name='inter_30')
	
	input_panda.to_csv(output_file+'_inter.tsv', sep='\t')
	input_panda_30.to_csv(output_file+'_inter_30.tsv', sep='\t')

if(datatype == "bulk_mega"):
	input_files = glob.glob(input_path + "/" + prefix + "*/aligned/inter.txt")
	input_files_30 = glob.glob(input_path + "/" + prefix + "*/aligned/inter_30.txt")
	
	input_data = {}
	input_data_30 = {}
	colnames = ["sample_name", "Sequenced Read Pairs", "Normal Paired", "Normal Paired (%)", "Chimeric Paired", "Chimeric Paired (%)", "Chimeric Ambiguous", "Chimeric Ambiguous (%)", "Unmapped", "Unmapped (%)", "Alignable (Normal+Chimeric Paired)", "Alignable (Normal+Chimeric Paired) (%)", "Unique Reads",  "PCR Duplicates", "Optical Duplicates", "Intra-fragment Reads", "Intra-fragment Reads (%/total)", "Intra-fragment Reads (%/unique)", "Below MAPQ Threshold", "Below MAPQ Threshold (%/total)", "Below MAPQ Threshold (%/unique)", "Hi-C Contacts", "Hi-C Contacts (%/total)", "Hi-C Contacts (%/unique)", "Ligation Motif Present", "Ligation Motif Present (%/total)", "Ligation Motif Present (%/unique)", "3' Bias (Long Range)", "Pair Type %(L-I-O-R)", "Inter-chromosomal", "Inter-chromosomal (%/total)", "Inter-chromosomal (%/unique)", "Intra-chromosomal", "Intra-chromosomal (%/total)", "Intra-chromosomal (%/unique)", "Short Range (<20Kb)", "Short Range (<20Kb) (%/total)", "Short Range (<20Kb) (%/unique)", "Long Range (>20Kb)", "Long Range (>20Kb) (%/total)", "Long Range (>20Kb) (%/unique)"]
	
	for input_file in input_files:
		sample_name = input_file.split('/')[-3]
		if sample_name == "mega":
			continue
		input_data[sample_name] = [sample_name]
		with open(input_file, 'r') as file:
			line_counter = 0
			for line in file:
				line_counter+=1
				entry = line.strip('\n').split(':')[1].split()
				if line_counter in [1, 7, 8, 9]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
				elif line_counter in [2, 3, 4, 5, 6]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
					input_data[sample_name].append(entry[1].split('(')[1].split(')')[0])
				elif line_counter in [10,11,12,13,16,17,18,19]:
					input_data[sample_name].append(int(entry[0].replace(',', '')))
					input_data[sample_name].append(entry[1].split('(')[1])
					input_data[sample_name].append(entry[3].split(')')[0])
				elif line_counter in [14, 15]:
					entry = line.strip('\n').split(':')[1]
					input_data[sample_name].append(entry)
	
	for input_file in input_files_30:
		sample_name = input_file.split('/')[-3]
		if sample_name == "mega":
			continue
		input_data_30[sample_name] = [sample_name]
		with open(input_file, 'r') as file:
			line_counter = 0
			for line in file:
				line_counter+=1
				entry = line.strip('\n').split(':')[1].split()
				if line_counter in [1, 7, 8, 9]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
				elif line_counter in [2, 3, 4, 5, 6]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
					input_data_30[sample_name].append(entry[1].split('(')[1].split(')')[0])
				elif line_counter in [10,11,12,13,16,17,18,19]:
					input_data_30[sample_name].append(int(entry[0].replace(',', '')))
					input_data_30[sample_name].append(entry[1].split('(')[1])
					input_data_30[sample_name].append(entry[3].split(')')[0])
				elif line_counter in [14, 15]:
					entry = line.strip('\n').split(':')[1]
					input_data_30[sample_name].append(entry)
	
	input_panda = pd.DataFrame.from_dict(input_data).T
	input_panda_30 = pd.DataFrame.from_dict(input_data_30).T
	
	input_panda.columns = colnames
	input_panda_30.columns = colnames
	
	input_panda = input_panda.sort_index()
	input_panda_30 = input_panda_30.sort_index()
	
	input_panda["Unique Reads (%/total)"] = input_panda["Unique Reads"]/input_panda["Sequenced Read Pairs"]*100
	input_panda["PCR Duplicates (%/total)"] = input_panda["PCR Duplicates"]/input_panda["Sequenced Read Pairs"]*100
	input_panda["Optical Duplicates (%/total)"] = input_panda["Optical Duplicates"]/input_panda["Sequenced Read Pairs"]*100
	
	input_panda_30["Unique Reads (%/total)"] = input_panda_30["Unique Reads"]/input_panda_30["Sequenced Read Pairs"]*100
	input_panda_30["PCR Duplicates (%/total)"] = input_panda_30["PCR Duplicates"]/input_panda_30["Sequenced Read Pairs"]*100
	input_panda_30["Optical Duplicates (%/total)"] = input_panda_30["Optical Duplicates"]/input_panda_30["Sequenced Read Pairs"]*100

	input_panda["Unique Reads (%/total)"] = input_panda["Unique Reads (%/total)"].astype(float).round(decimals=2)
	input_panda["PCR Duplicates (%/total)"] = input_panda["PCR Duplicates (%/total)"].astype(float).round(decimals=2)
	input_panda["Optical Duplicates (%/total)"] = input_panda["Optical Duplicates (%/total)"].astype(float).round(decimals=2)
	
	input_panda_30["Unique Reads (%/total)"] = input_panda_30["Unique Reads (%/total)"].astype(float).round(decimals=2)
	input_panda_30["PCR Duplicates (%/total)"] = input_panda_30["PCR Duplicates (%/total)"].astype(float).round(decimals=2)
	input_panda_30["Optical Duplicates (%/total)"] = input_panda_30["Optical Duplicates (%/total)"].astype(float).round(decimals=2)

	input_panda["Unique Reads (%/total)"] = input_panda["Unique Reads (%/total)"].astype(str) + "%"
	input_panda["PCR Duplicates (%/total)"] = input_panda["PCR Duplicates (%/total)"].astype(str) + "%"
	input_panda["Optical Duplicates (%/total)"] = input_panda["Optical Duplicates (%/total)"].astype(str) + "%"

	input_panda_30["Unique Reads (%/total)"] = input_panda_30["Unique Reads (%/total)"].astype(str) + "%"
	input_panda_30["PCR Duplicates (%/total)"] = input_panda_30["PCR Duplicates (%/total)"].astype(str) + "%"
	input_panda_30["Optical Duplicates (%/total)"] = input_panda_30["Optical Duplicates (%/total)"].astype(str) + "%"
	
	input_panda = input_panda[["sample_name", "Sequenced Read Pairs", "Normal Paired", "Normal Paired (%)", "Chimeric Paired", "Chimeric Paired (%)", "Chimeric Ambiguous", "Chimeric Ambiguous (%)", "Unmapped", "Unmapped (%)", "Alignable (Normal+Chimeric Paired)", "Alignable (Normal+Chimeric Paired) (%)", "Unique Reads", "Unique Reads (%/total)",  "PCR Duplicates", "PCR Duplicates (%/total)", "Optical Duplicates", "Optical Duplicates (%/total)", "Intra-fragment Reads", "Intra-fragment Reads (%/total)", "Intra-fragment Reads (%/unique)", "Below MAPQ Threshold", "Below MAPQ Threshold (%/total)", "Below MAPQ Threshold (%/unique)", "Hi-C Contacts", "Hi-C Contacts (%/total)", "Hi-C Contacts (%/unique)", "Ligation Motif Present", "Ligation Motif Present (%/total)", "Ligation Motif Present (%/unique)", "3' Bias (Long Range)", "Pair Type %(L-I-O-R)", "Inter-chromosomal", "Inter-chromosomal (%/total)", "Inter-chromosomal (%/unique)", "Intra-chromosomal", "Intra-chromosomal (%/total)", "Intra-chromosomal (%/unique)", "Short Range (<20Kb)", "Short Range (<20Kb) (%/total)", "Short Range (<20Kb) (%/unique)", "Long Range (>20Kb)", "Long Range (>20Kb) (%/total)", "Long Range (>20Kb) (%/unique)"]]
	input_panda_30 = input_panda_30[["sample_name", "Sequenced Read Pairs", "Normal Paired", "Normal Paired (%)", "Chimeric Paired", "Chimeric Paired (%)", "Chimeric Ambiguous", "Chimeric Ambiguous (%)", "Unmapped", "Unmapped (%)", "Alignable (Normal+Chimeric Paired)", "Alignable (Normal+Chimeric Paired) (%)", "Unique Reads", "Unique Reads (%/total)",  "PCR Duplicates", "PCR Duplicates (%/total)", "Optical Duplicates", "Optical Duplicates (%/total)", "Intra-fragment Reads", "Intra-fragment Reads (%/total)", "Intra-fragment Reads (%/unique)", "Below MAPQ Threshold", "Below MAPQ Threshold (%/total)", "Below MAPQ Threshold (%/unique)", "Hi-C Contacts", "Hi-C Contacts (%/total)", "Hi-C Contacts (%/unique)", "Ligation Motif Present", "Ligation Motif Present (%/total)", "Ligation Motif Present (%/unique)", "3' Bias (Long Range)", "Pair Type %(L-I-O-R)", "Inter-chromosomal", "Inter-chromosomal (%/total)", "Inter-chromosomal (%/unique)", "Intra-chromosomal", "Intra-chromosomal (%/total)", "Intra-chromosomal (%/unique)", "Short Range (<20Kb)", "Short Range (<20Kb) (%/total)", "Short Range (<20Kb) (%/unique)", "Long Range (>20Kb)", "Long Range (>20Kb) (%/total)", "Long Range (>20Kb) (%/unique)"]]
	
	with pd.ExcelWriter(output_file+'_qc.xlsx') as writer:  
	    input_panda.to_excel(writer, sheet_name='inter')
	    input_panda_30.to_excel(writer, sheet_name='inter_30')
	
	input_panda.to_csv(output_file+'_inter.tsv', sep='\t')
	input_panda_30.to_csv(output_file+'_inter_30.tsv', sep='\t')





