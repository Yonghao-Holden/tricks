#!/usr/bin/env python3
"""
Make session file and upload for washu browser view

Usage: python3 /bar/yliang/tricks/makesession.py --input ${INPUT} --output ${OUTPUT} --prefix ${PREFIX} --suffix ${SUFFIX} --bundle ${BUNDLE}
"""

import argparse
import sys
import os
import glob


# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "path to input files", type = str, default = "/bar/yliang/public_html/Epitherapy_3D/rna-seq")
parser.add_argument("--output", help = "output file basename", type = str, default = "browser_session")
parser.add_argument("--prefix", help = "prefix of input files", type = str, default = "")
parser.add_argument("--suffix", help = "suffix of input files", type = str, default = "")
parser.add_argument("--bundle", help = "bundle number", type = str, default = "aa0d66c0-cfb4-11eb-85ca-67b9cb6fe783")

args = parser.parse_args()

input_path = args.input
output_file = args.output + ".json"
prefix = args.prefix
suffix = args.suffix
bundle = args.bundle

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]

def remove_suffix(text, suffix):
	return text.strip(suffix)

###############################################################
# Part 1. Data Read in
###############################################################
input_files = glob.glob(input_path + "/" + prefix + "*" + suffix)

input_files = [x for x in input_files if x[-3:] != "bdg" and x[-8:] != "bedgraph"]

input_files = sorted(input_files, key=str.lower)

#input_file_names = [remove_prefix(remove_suffix(os.path.basename(file), suffix), prefix) for file in input_files]
input_file_names = [os.path.basename(file) for file in input_files]

input_files = ["https://wangftp.wustl.edu/~yliang/" + remove_prefix(x, "/bar/yliang/public_html/") for x in input_files]

output = []

output.append("""{"genomeName":"hg38","viewInterval":{"start":1259057700,"end":1259378068},"tracks":[{"name":"refGene","type":"geneannotation","label":"refGene","options":{"label":"refGene"},"url":"","metadata":{"Track type":"geneannotation"},"isSelected":false,"fileObj":"","files":[],"tracks":[],"isText":false,"textConfig":{},"genome":"hg38"},{"name":"gencodeV29","type":"geneannotation","label":"gencodeV29","options":{"label":"gencodeV29"},"url":"","metadata":{"Track type":"geneannotation"},"isSelected":false,"fileObj":"","files":[],"tracks":[],"isText":false,"textConfig":{},"genome":"hg38"},{"name":"RepeatMasker","type":"repeatmasker","label":"RepeatMasker","options":{"label":"RepeatMasker"},"url":"https://vizhub.wustl.edu/public/hg38/rmsk16.bb","metadata":{"Track type":"repeatmasker"},"isSelected":false,"fileObj":"","files":[],"tracks":[],"isText":false,"textConfig":{}},{"name":"Ruler","type":"ruler","label":"Ruler","options":{"label":"Ruler"},"url":"","metadata":{"Track type":"ruler"},"isSelected":false,"fileObj":"","files":[],"tracks":[],"isText":false,"textConfig":{}}""")
for name, file in zip(input_file_names, input_files):
	if "tbi" in file:
		continue
	elif "bdg" in file or "bedgraph" in file:
		file_type = "bedgraph"
	elif "bed" in file or "narrowPeak" in file:
		file_type = "bed"
	elif "bigwig" in file or "bw" in file:
		file_type = "bigwig"
	output.append(""",{"name":"%s","type":"%s","label":"%s","options":{"label":"%s"},"url":"%s","metadata":{"Track type":"%s"},"isSelected":false,"fileObj":"","files":[],"tracks":[],"isText":false,"textConfig":{},"datahub":"Custom track","selectedTabIndex":0,"trackAdded":false,"urlError":""}""" %(name, file_type, name, name, file, file_type))
output.append("""],"metadataTerms":[],"regionSets":[],"regionSetViewIndex":-1,"trackLegendWidth":120,"bundleId":"%s","isShowingNavigator":true,"isShowingVR":false,"layout":{}}""" %(bundle))


with open(output_file, 'w+') as file:
	file.write(''.join(output))


