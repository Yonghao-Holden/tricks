#!/usr/bin/env python3
"""
Input mega data csv file (Label_to_Sample_Mapping_File in the meta documentation tab) downloaded from CPTAC and generate annotation.txt files in each folder for philosopher pipeline

Need to have mannual curation on the csv file first before putting it through this script. 
Check (1) folder name, (2) check reference sample, (3) check disqualified sample, (4) remove comments at the bottom

Usage: python3 CPTAC_megadata_parser.py --input <mega.csv>


"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import re
import warnings
import math
warnings.simplefilter("ignore")

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input file", type = str)
parser.add_argument("--data", help = "data type", type = str)
parser.add_argument("--exp", help = "number of experiment", type = int, default=12)
args = parser.parse_args()

input_file = args.input
data_type = args.data
number_of_exp = args.exp

###############################################################
# Part 1. Data Read in
############################################################### 
def extract_information(input_row, data_type):
	folder_name = str(input_row["Folder Name"])
	if data_type == "tmt11":
		useful_columns = ["TMT11-126 Specimen Label","TMT11-127N Specimen Label","TMT11-127C Specimen Label","TMT11-128N Specimen Label","TMT11-128C Specimen Label",\
		"TMT11-129N Specimen Label","TMT11-129C Specimen Label","TMT11-130N Specimen Label","TMT11-130C Specimen Label","TMT11-131N Specimen Label","TMT11-131C Specimen Label"]
		labels = ["126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C"]
	if data_type == "tmt10":
		useful_columns = ["TMT10-126 Specimen Label","TMT10-127N Specimen Label","TMT10-127C Specimen Label","TMT10-128N Specimen Label","TMT10-128C Specimen Label",\
		"TMT10-129N Specimen Label","TMT10-129C Specimen Label","TMT10-130N Specimen Label","TMT10-130C Specimen Label","TMT10-131 Specimen Label"]
		labels = ["126","127N","127C","128N","128C","129N","129C","130N","130C","131N"]

	with open(folder_name + "/annotation.txt", 'w') as file:
		for useful_column, label in zip(useful_columns, labels):
			#print(str(input_row[useful_column]))
			if str(input_row[useful_column]) != "nan":
				file.write(label+" "+str(input_row[useful_column])+"\n")
			elif str(input_row[useful_column]) == "nan":
				file.write(label+" "+str(input_row[useful_column.replace("Specimen Label","Participant ID")])+"\n")

input_panda = pd.read_excel(input_file)
input_panda = input_panda.head(number_of_exp)
input_panda.apply(extract_information, data_type=data_type,axis=1)






