#!/usr/bin/env python3
"""
Run it in a CPTAC mass spec data folder and will genearte a manifest file for fragpipe run

path needs to be absolute path

Usage: python3 CPTAC_megadata_parser.py --study <study name>
"""

import argparse
import sys
import os
import numpy as np
import glob

# Get the path to input files
parser = argparse.ArgumentParser()
#parser.add_argument("--input", help = "input path name", type = str)
parser.add_argument("--study", help = "study name", type = str)
parser.add_argument("--exptype", help = "experiment type for MS, DDA or DIA, check on CPTAC website, Protocol tab", type = str, default = "DDA")
args = parser.parse_args()

#input_path = args.input
input_path = os.getcwd()
study = args.study
exp_type = args.exptype

###############################################################
# Part 1. Data Read in
###############################################################
input_files = glob.glob(input_path+"/*/*mzML")
output_string = ""

for input_file in input_files:
    folder_name = input_file.split("/")[-2]
    output_string += input_file +'\t'+ folder_name + "\t\t" + exp_type + "\n"

with open(study + "_manifest.fp-manifest", 'w') as file:
    file.write(output_string)
 