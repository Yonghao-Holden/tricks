#!/bin/bash
module load juicer
#SBATCH --time=10-00:00:00

#find ./*/aligned -name merged_nodups.txt.gz | while read file; do echo "gunzip $file" >> gunzip_commands.txt ; done;
#JOB=`ll ./*/aligned/merged_nodups.txt.gz | wc -l`
#parallel_GNU -j $JOB < gunzip_commands.txt

#/scratch/twlab/mayank-choudhary/jobs/htcf_mega10.sh -g hg38 -s MboI -d $1
/home/yonghao/softwares/juicer/SLURM/scripts/mega_10days.sh -g hg38 -s MboI -d $1

#find ./*/aligned -name merged_nodups.txt | while read file; do echo "gzip $file" >> gzip_commands.txt ; done;
#JOB=`ll ./*/aligned/merged_nodups.txt | wc -l`
#parallel_GNU -j $JOB < gzip_commands.txt