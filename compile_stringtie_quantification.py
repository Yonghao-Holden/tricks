#!/usr/bin/env python3
##########################################################################################################
# Author: Holden Liang
# this script is to compile gene expression file downloaded from GDC together to one dataframe (tsv)
##########################################################################################################
import argparse
import glob
import pandas as pd
from tqdm import tqdm

#########################################
## 1. parser
#########################################
parser = argparse.ArgumentParser()
parser.add_argument("--input", help="path to gene expression files",type = str, default="./*stats/t_data.ctab")
parser.add_argument("--output", help="output file name",type = str, default="rna-seq_dbgap_phs000892_CPTAC2_stringtie_TSTET.txt")
args = parser.parse_args()

input_path = args.input
outpu_file_name = args.output

input_files = glob.glob(input_path)

output_panda = pd.DataFrame(columns=["t_name","gene_name","FPKM","sample"])

for input_file in tqdm(input_files):
	input_panda = pd.read_table(input_file)
	## ./ffc4018e-df78-4b25-9de4-91ec58ca5b1f.rna_seq.genomic.gdc_realn.stats/t_data.ctab
	input_panda["sample"] = input_file.split("/")[1].replace("_R1.fastqAligned.sortedByCoord.out.stats","").replace("Trimmed_","")
	output_panda = pd.concat([output_panda, input_panda[["t_name","gene_name","FPKM","sample"]]])

output_panda = output_panda.pivot(index=["t_name","gene_name"], columns="sample")
output_panda.columns = output_panda.columns.droplevel(0)
output_panda.columns.name = None

output_panda.to_csv(outpu_file_name, sep="\t")



