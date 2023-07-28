#!/bin/bash
# programmer: Holden Liang
# Run it like this: bash /home/yonghao/tricks/wgbs_pipe_htcf_starter_v1.sh -i <input dir> 
# Run it in the base folder (~/wgbs), under which is ~/wgbs/fastq. Put ALL fastq files under this fastq folder.
# it will run all samples in the folder at the same time. (submit jobs for all samples together, so don't put too many samples in the folder)
###################################################
# run this script on htcf to run bismark pipeline, will generate methylation level for CpG and basic qc plots including correlation heatmap across samples and methylation level distribution.
# Written by Yonghao Liang (Holden)
# for sbatch flags and how to run jobs one by one, check out websites below:
# https://slurm.schedmd.com/sbatch.html
# https://hpc.nih.gov/docs/job_dependencies.html
###################################################

# bash /home/yonghao/tricks/wgbs_pipe_htcf_starter_v4.sh -i . -a yes

#########################################################################################################################################################
## INPUT
# Arguments to bash script to decide on various parameters
INPUT_DIR="." # it should have a fastq folder in it
READ1EXTENSION="_R1.fastq.gz"
GENOME="hg38"
MAXJOBSALIGN=12 # Cores for each alignment
## only switch $ALLTOGETHER to yes if you have less then 3 samples to run. 
## If yes, it will submit jobs for ALL samples in your fastq folder to htcf, and each sample will have multiple jobs submitted at the same time.
## so it can easily take over the whole server if you have too many samples
ALLTOGETHER=yes 
PBAT=no
#########################################################################################################################################################

while [[ $# > 1 ]]
do
key="$1"
case $key in
    -i|--input)
    INPUT_DIR="$2"
    shift # past argument
    ;;
    -a|--alltogether)
    ALLTOGETHER="$2"
    shift # past argument
    ;;
    -p|--pbat)
    PBAT="$2"
    shift # past argument
    ;;
    -r|--read1extension)
    READ1EXTENSION="$2"
    shift # past argument
    ;;
    -g|--genome)
    GENOME="$2"
    shift # past argument
    ;;
    -m|--max)
    MAXJOBSALIGN="$2"
    shift # past argument
    ;;    
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

## on htcf(SLURM)
# squeue -u yonghao
##############################################################
# 1. copy fastq files and genome annotation files from our sever to htcf
# rsync -avFP /path/to/file yonghao@65.254.100.82:/path/to/file
# parallel_GNU -j 10 < fastq_rsync.sh

# rsync -avFPL /bar/yliang/genomes/public/hg38/bismark yonghao@65.254.100.82:/scratch/twlab/yliang/genome/hg38
# rsync -avFPL /bar/yliang/genomes/public/mm10/bismark yonghao@65.254.100.82:/scratch/twlab/yliang/genome/mm10
# rsync -avFPL /bar/yliang/genomes/public/lambda yonghao@65.254.100.82:/scratch/twlab/yliang/genome

# find /scratch/twlab/yliang/genome/* | while read file; do touch $file; done

##############################################################
# 2. combine 2 fastq files generated in 2 different lanes (optional)
# zcat sample.Lane1.R1.fastq.gz sample.Lane2.R1.fastq.gz | gzip > sample.R1.fastq.gz
# parallel_GNU -j 10 < fastq_merge.sh
# rm *Lane*

# Example:
# rsync -avFP yliang@10.20.127.3:/taproom/data/dli/NovaSeq-S4-DT-2861687/xfer.genome.wustl.edu/gxfer1/20072858207592/CTTGGAA_S16_L003_R1_001.fastq.gz HiC_93VU147T_HPV_BRep1_TRep1_Lane1_R1.fastq.gz 
# zcat HiC_93VU147T_HPV_BRep1_TRep1_Lane1_R1.fastq.gz HiC_93VU147T_HPV_BRep1_TRep1_Lane2_R1.fastq.gz | gzip > HiC_93VU147T_HPV_BRep1_TRep1_R1.fastq.gz

##############################################################
# 3. run bismark
# (1) Get scripts from our sever
# rsync -avFP /bar/yliang/tricks/wgbs_pipe_htcf*.sh yonghao@65.254.100.82:/home/yonghao/tricks/
# bash /home/yonghao/tricks/wgbs_pipe_htcf_starter_v3.sh -i . -a yes

# (2) Run bismark on each replicate and clean up (submit job one by one based on the generation of fincin.out)
cd $INPUT_DIR
if [[ $ALLTOGETHER != "yes" ]]; then
    for FILE1 in ./fastq/*${READ1EXTENSION}
    do
        FILE2=${FILE1/R1/R2}
        xbase="$(basename $FILE1 $READ1EXTENSION)"
        echo $xbase
        mkdir $xbase $xbase/fastq
        cd $xbase/fastq
        ln -s ../../$FILE1 .
        ln -s ../../$FILE2 .
        cd ..
        echo "start $xbase bismark run"
        bash /home/yonghao/tricks/wgbs_pipe_htcf_v4.1.sh --sample $xbase --genome $GENOME --input $PWD --max $MAXJOBSALIGN --pbat $PBAT
        while true; do
            #FILENAME=`inotifywait -m --quiet --recursive -e create,modify,close,move --format '%f' ./`
            #if [[ $FILENAME =~ done-[0-9]{8}.out ]]; then
            if [[ `ls ./debug/done-*.out > /dev/null 2> /dev/null; echo $?` -eq 0 ]]; then
                echo "$xbase bismark run finished"
                break
            fi
        done
        cd ..
    done
else
    for FILE1 in ./fastq/*${READ1EXTENSION}
    do
        FILE2=${FILE1/R1/R2}
        xbase="$(basename $FILE1 $READ1EXTENSION)"
        echo $xbase
        mkdir $xbase $xbase/fastq
        cd $xbase/fastq
        ln -s ../../$FILE1 .
        ln -s ../../$FILE2 .
        cd ..
        echo "start $xbase bismark run"
        bash /home/yonghao/tricks/wgbs_pipe_htcf_v4.1.sh --sample $xbase --genome $GENOME --input $PWD --max $MAXJOBSALIGN --pbat $PBAT
        # sbatch /scratch/twlab/yliang/juicer_scripts/run-juicer-hg38.sbatch HiC_93VU147T_HPV_BRep1_TRep1
        cd ..
    done    
fi
## when the pipeline is finished with $ALLTOGETHER=="no"
## you should have cov.gz files, correlation heatmap, methylation level distribution plots in your folder
## and run multiqc .
#if [[ $ALLTOGETHER == "no" ]]; then
#    multiQC .
#fi




