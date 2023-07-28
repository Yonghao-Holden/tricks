#!/usr/bin/env python3
"""
Input fasta file of cDNA, output fasta file of all ORFs that are longer than <threshold> amino acid

Usage: python3 cDNAtoprotein.py --input <cDNA.fa> --length 20 --output <protein.fa>
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import re

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "path to input cDNA fasta file", type = str)
parser.add_argument("--length", help = "length of amino acid", type = int)
parser.add_argument("--output", help = "name of output protein sequence fasta file", type = str)
args = parser.parse_args()

input_file = args.input
length = args.length
output_file = args.output


###############################################################
# Part 1. Functions
###############################################################
def readfile(input_file):
    output_dict = {}
    with open(input_file, 'r') as file:
        for line in file:
            entry = line.strip('\n')
            if entry[0] == ">":
                output_dict[entry.split(">")[1]] = ""
                transcript_name = entry.split(">")[1]
            else:
                output_dict[transcript_name]+=entry
    return(output_dict)
    
def translate(seq):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'!', # ! represent start codon
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    for i in range(0, len(seq), 3): # len(seq)=10. 0, 3, 6, 9
        if i+3 <= len(seq):
            codon = seq[i:i + 3]
            protein+= codon_table[codon]
    return protein

#def find_start_codon(seq):
#    location = {"1"=[],"2"=[],"3"=[]}
#    for i in range(len(seq)):
#        if seq[i:i+3] == "ATG":
#            location[str(i/3)].append(i)
#
#def find_stop_codon(seq):
#    location = {"1"=[],"2"=[],"3"=[]}
#    for i in range(len(seq)):
#        if seq[i:i+3] == "TAG" or seq[i:i+3] == "TAA" or seq[i:i+3] == "TGA":
#            location[str(i/3)].append(i)

###############################################################
# Part 2. Main
###############################################################
input_data = readfile(input_file)
#print(input_data)
output_data = {}
transcript_number = 1

for transcript in input_data:
    cDNA_sequence = input_data[transcript]
    #print(transcript_number, end =" ")
    transcript_number += 1
    sequence_number = 1
    for k in [0,1,2]: # three ORFs
        protein_sequence = translate(cDNA_sequence[k:])
        #print(k)
        #print("cDNA sequence")
        #print(cDNA_sequence[k:])
        #print("protein sequence")
        #print(protein_sequence)
        for i in range(len(protein_sequence)): # len(protein_sequence) = 10 (0,1,2...9)
            if protein_sequence[i] == "!": # i = 3
                for j in range(i, len(protein_sequence)): # j in range(3,10), ie 3,4...9
                    if protein_sequence[j] == "_":  # j = 7
                        if len(protein_sequence[i:j]) < length:
                            break
                        else:
                            output_data[transcript+"_protein_sequence_"+str(sequence_number)] = protein_sequence[i:j].replace("!", "M") # [3,7]
                            sequence_number+=1
                            break
                    elif j == len(protein_sequence)-1: # j = 9
                        if len(protein_sequence[i:j+1]) < length:
                            break
                        else:
                            output_data[transcript+"_protein_sequence_"+str(sequence_number)] = protein_sequence[i:j+1].replace("!", "M") # M + [3,10]
                            sequence_number+=1
                            break

output_string = ""
for transcript in output_data:
    output_string += ">tr|" + transcript + "|" + transcript + " " + transcript + " " + "OS=Homo sapiens OX=9606 GN=" + transcript.split("_")[2] +"\n"+output_data[transcript]+"\n"

with open(output_file, 'w') as file:
    file.write(output_string)

