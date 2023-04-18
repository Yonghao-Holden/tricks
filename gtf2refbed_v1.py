# TEProf3: TE-derived Promoter Finder 3
# Author: Holden Liang
# visualizing gtf files for epigenome browser
# convert gtf files from stringtie to refbed file
#------------------------------------------------------------------------
import argparse
import subprocess # https://geekflare.com/python-run-bash/
import os
import multiprocess as mp

parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input file", type = str)

args = parser.parse_args()

input_file = args.input

def gtf_to_refbed(input_gtf_file):
	""" convert stringtie gtf file to refbed file for washu genome browser visualization """
	# chr, transcript_start, transcript_stop, translation_start, translation_stop, 
	# strand, gene_name, transcript_id, type, exon(including UTR bases) starts, 
	# exon(including UTR bases) stops, and additional gene info (optional)
	output_refbed_file = input_gtf_file.replace("gtf", "refbed")
	refbed_data = {}

	with open(input_gtf_file, 'r') as input_file:
		last_line = input_file.readlines()[-1]

	with open(input_gtf_file, 'r') as input_file:
		input_file.readline()
		input_file.readline()
		first_transcript = 1
		for line in input_file:
			entry = line.strip('\n').split('\t')
			if entry[2] == "transcript":
				if first_transcript == 1:
					first_transcript = 0
				elif first_transcript == 0:
					refbed_data[transcript_id] += [','.join(exon_starts), ','.join(exon_stops) ,detail_info]
				detail_info = entry[8]
				gene_id = detail_info.split("gene_id \"")[1].split("\";")[0]
				transcript_id = detail_info.split("transcript_id \"")[1].split("\";")[0]
				refbed_data[transcript_id] = [entry[0], entry[3], entry[4], entry[3], entry[4], entry[6], gene_id, transcript_id, "coding"]
				exon_starts = []
				exon_stops = []
			elif entry[2] == "exon":
				exon_starts.append(entry[3])
				exon_stops.append(entry[4])
				if last_line is line:
					refbed_data[transcript_id] += [','.join(exon_starts), ','.join(exon_stops) ,detail_info]
	with open(output_refbed_file, 'w') as output_file:	
		for transcript_id in refbed_data:
			output_file.write("\t".join(refbed_data[transcript_id])+'\n')
	subprocess.run("sort -k1,1 -k2,2n "+ str(output_refbed_file) +" | bgzip > "+ str(output_refbed_file) +".sorted.gz", shell=True)
	subprocess.run("tabix -p bed "+ str(output_refbed_file) +".sorted.gz", shell=True)

gtf_to_refbed(input_file)



