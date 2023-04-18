#!/usr/bin/env python3
"""
Input mega data csv file downloaded from CPTAC and generate annotation.txt files in each folder for philosopher pipeline

Usage: python3 CPTAC_megadata_parser.py --input <mega.csv>
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import re

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input path", type = str)
args = parser.parse_args()

input_file = args.input

###############################################################
# Part 1. Data Read in
###############################################################
def readfile(input_file):
    with open(input_file, 'r') as file:
        header = file.readline().strip('\n').split(',')
        header = header[3:]
        allinput = ""
        for line in file:
            entry = line.strip('\n')
            allinput += entry + ' '
    allinput = re.split('","|""|"', allinput)
    #print(allinput)
    allinput = remove_space(allinput)
    #print(allinput)
    return(header, allinput)
    
def remove_space(input_list):
    output_list = []
    for i in input_list:
        if str(i) == "" or str(i) == " ":
            continue
        i = str(i).replace(' ','_')
        output_list.append(i)
    return(output_list)

# https://github.com/Nesvilab/philosopher/blob/47ad5f7394ad8aa563c0396ea2f9f144bdf7af20/lib/qua/qua.go
def convert_channel_name(channel_name):
    if channel_name == "tmt_126":
        return("126")
    elif channel_name == "tmt_127n":
        return("127N")
    elif channel_name == "tmt_127c":
        return("127C")
    elif channel_name == "tmt_128n":
        return("128N")
    elif channel_name == "tmt_128c":
        return("128C")
    elif channel_name == "tmt_129n":
        return("129N")
    elif channel_name == "tmt_129c":
        return("129C")
    elif channel_name == "tmt_130n":
        return("130N")
    elif channel_name == "tmt_130c":
        return("130C")
    elif channel_name == "tmt_131":
        return("131N")
    elif channel_name == "tmt_131C":
        return("131C")
    elif channel_name == "itraq_114":
        return("114")
    elif channel_name == "itraq_115":
        return("115")
    elif channel_name == "itraq_116":
        return("116")
    elif channel_name == "itraq_117":
        return("117")       
    
header, input_data = readfile(input_file)
#print(header, input_data)
    
transformed_input_data = {}

for i in range(len(input_data)):
    col_number = i % (len(header)+3)
    if col_number == 1:
        dataset = '_'.join(re.split(' ', input_data[i]))
        transformed_input_data[dataset] = []
        #print(dataset)
    elif col_number == 0 or col_number == 2:
        continue
    else:
        transformed_input_data[dataset].append(input_data[i])
        
#print(transformed_input_data)
    
for dataset in transformed_input_data:
    output_string = ""
    for channel_name, sample in zip(header, transformed_input_data[dataset]):
        channel_name_converted = convert_channel_name(channel_name)
        output_string += str(channel_name_converted) +' '+ str(sample) + "\n"
    with open(dataset + "/annotation.txt", 'a+') as file:
        file.write(output_string)
    print(output_string)
       


