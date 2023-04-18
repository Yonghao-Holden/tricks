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

def get_simplified_transcript_type(transcript_type):
	transcript_type_list_coding = ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "polymorphic_pseudogene", "protein_coding", "IG_LV_gene"]
	transcript_type_list_pseudo = ["rRNA_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "pseudogene", "processed_pseudogene", "processed_pseudogene"]
	transcript_type_list_nonCoding = ["non_coding", "antisense_RNA", "ribozyme", "sRNA", "3prime_overlapping_ncRNA", "macro_lncRNA", "bidirectional_promoter_lncRNA", "scaRNA", "scRNA", "snoRNA", "snRNA", "sense_overlapping", "sense_intronic", "3prime_overlapping_ncrna", "Mt_rRNA", "Mt_tRNA", "antisense", "lincRNA", "miRNA", "misc_RNA", "processed_transcript","rRNA","lncRNA"]
	transcript_type_list_problem = ["TEC"]
	if transcript_type in transcript_type_list_coding:
		return("coding")
	elif transcript_type in transcript_type_list_pseudo:
		return("pseudo")
	elif transcript_type in transcript_type_list_nonCoding:
		return("nonCoding")
	elif transcript_type in transcript_type_list_problem:
		return("problem")
	else:
		return("other")

def gtf_to_refbed(input_gtf_file):
	""" convert stringtie gtf file to refbed file for washu genome browser visualization """
	# chr, transcript_start, transcript_stop, translation_start, translation_stop, 
	# strand, gene_name, transcript_id, type, exon(including UTR bases) starts, 
	# exon(including UTR bases) stops, and additional gene info (optional)

	subprocess.run("awk '$3~/transcript/ {print $0}' "+input_gtf_file+" > "+input_gtf_file.replace("gtf","temp.gtf"), shell=True)
	subprocess.run("awk '($3~/exon/ || $3~/start_codon/ || $3~/stop_codon/) {print $0}' "+input_gtf_file+" >> "+input_gtf_file.replace("gtf","temp.gtf"), shell=True)

	output_refbed_file = input_gtf_file.replace("gtf", "refbed")
	input_data = {}
	refbed_data = {}

	useful_info = ["transcript", "exon", "start_codon", "stop_codon"]

	with open(input_gtf_file.replace("gtf","temp.gtf"), 'r') as input_file:
		for line in input_file:
			entry = line.strip('\n').split('\t')
			detail_info = entry[8]
			if entry[2] in useful_info:
				transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
				try:
					gene_name = detail_info.split("gene_id \"")[1].split("\";")[0]
				except:
					gene_name = transcript_id
				try:
					gene_name = detail_info.split("gene_name \"")[1].split("\";")[0]
				except:
					gene_name = transcript_id
				try:
					transcript_type = get_simplified_transcript_type(detail_info.split("transcript_biotype \"")[1].split("\";")[0])
				except:
					transcript_type = "coding"
				try:
					transcript_type = get_simplified_transcript_type(detail_info.split("gene_type \"")[1].split("\";")[0])
				except:
					transcript_type = "coding"
			if entry[2] == "transcript":
				input_data[transcript_id] = [[entry[0], entry[3], entry[4], entry[3], entry[4], entry[6], gene_name, transcript_id, transcript_type],[],[],[detail_info]]
			elif entry[2] == "exon":
				input_data[transcript_id][1].append(entry[3])
				input_data[transcript_id][2].append(entry[4])
			elif entry[2] == "start_codon":
				input_data[transcript_id][0][3] = entry[3]
			elif entry[2] == "stop_codon":
				input_data[transcript_id][0][4] = str(int(entry[3])-1)

	with open(output_refbed_file, 'w') as output_file:	
		for transcript_id in input_data:
			output_file.write('\t'.join(input_data[transcript_id][0] + [','.join(input_data[transcript_id][1])] + [','.join(input_data[transcript_id][2])] + input_data[transcript_id][3])+'\n')

	subprocess.run("sort -k1,1 -k2,2n "+ str(output_refbed_file) +" | bgzip > "+ str(output_refbed_file) +".sorted.gz", shell=True)
	subprocess.run("tabix -p bed "+ str(output_refbed_file) +".sorted.gz", shell=True)
	subprocess.run("rm "+input_gtf_file.replace("gtf","temp.gtf"), shell=True)


gtf_to_refbed(input_file)



