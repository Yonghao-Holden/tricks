#!/usr/bin/env python3
'''
# Author: Holden Liang
# visualizing gtf files for epigenome browser
# convert gtf files from stringtie, TACO to refbed file (only have transcript and exon at column 3)
# run it as followed:
# /bar/yliang/tricks/gtf2refbed_v2.py --input <input.gtf>
'''
import argparse
import subprocess # https://geekflare.com/python-run-bash/
import os
import multiprocess as mp

parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input file", type = str, default="yo")

args = parser.parse_args()

input_file = args.input

if input_file == "yo":
	print("please include input.gtf file. run it as followed:")
	print("/bar/yliang/tricks/gtf2refbed_v2.py --input <input.gtf>")
	exit()

def gtf_to_refbed(input_gtf_file):
	""" convert stringtie gtf file to refbed file for washu genome browser visualization """
	# chr, transcript_start, transcript_stop, translation_start, translation_stop, 
	# strand, gene_name, transcript_id, type, exon(including UTR bases) starts, 
	# exon(including UTR bases) stops, and additional gene info (optional)
	output_refbed_file = input_gtf_file.replace("gtf", "refbed")
	input_data = {}
	refbed_data = {}

	with open(input_gtf_file, 'r') as input_file:
		for line in input_file:
			if line[0] != "#":
				entry = line.strip('\n').split('\t')
				detail_info = entry[8]
				if entry[2] == "transcript" or entry[2] == "exon":
					gene_id = detail_info.split("gene_id \"")[1].split("\";")[0]
					transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
				if entry[2] == "transcript":
					if transcript_id in input_data:
						input_data[transcript_id][0] = [entry[0], entry[3], entry[4], entry[3], entry[4], entry[6], gene_id, transcript_id, "coding"]
						input_data[transcript_id][3] = [detail_info]
					if transcript_id not in input_data:
						input_data[transcript_id] = [[entry[0], entry[3], entry[4], entry[3], entry[4], entry[6], gene_id, transcript_id, "coding"],[],[],[detail_info]]
				elif entry[2] == "exon":
					if transcript_id not in input_data:
						input_data[transcript_id] = [[],[],[],[]]
					input_data[transcript_id][1].append(entry[3])
					input_data[transcript_id][2].append(entry[4])

	with open(output_refbed_file, 'w') as output_file:	
		for transcript_id in input_data:
			output_file.write('\t'.join(input_data[transcript_id][0] + [','.join(input_data[transcript_id][1])] + [','.join(input_data[transcript_id][2])] + input_data[transcript_id][3])+'\n')

	subprocess.run("sort -k1,1 -k2,2n "+ str(output_refbed_file) +" | bgzip > "+ str(output_refbed_file) +".sorted.gz", shell=True)
	subprocess.run("tabix -p bed "+ str(output_refbed_file) +".sorted.gz", shell=True)

gtf_to_refbed(input_file)



