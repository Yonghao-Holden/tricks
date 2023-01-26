#!/bin/bash
#!/usr/bin/awk
# programmer: Holden Liang
# Run it in the base folder (~/Hi-C), under which is ~/Hi-C/fastq. Put ALL fastq files under this fastq folder.
# bash /home/yonghao/tricks/juicer_htcf_v2.sh -i /scratch/twlab/yliang/HNSCC/data/hic/UMSCC47 -r _Lane1_R1.fastq.gz -g hg38 -a yes
###################################################
# run this script on htcf to generate .hic files using juicer
# Written by Yonghao Liang(Holden)
# for sbatch flags and how to run jobs one by one, check out websites below:
# https://slurm.schedmd.com/sbatch.html
# https://hpc.nih.gov/docs/job_dependencies.html
###################################################

#########################################################################################################################################################
## INPUT
# Arguments to bash script to decide on various parameters
INPUT_DIR="." # it should have a fastq folder in it
READ1EXTENSION="_R1.fastq.gz"
GENOME="hg38"
ENZYME="MboI"

ALLTOGETHER=no 

#########################################################################################################################################################

while [[ $# > 1 ]]
do
key="$1"
case $key in
    -i|--input)
    INPUT_DIR="$2"
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
    -e|--enzyme)
    ENZYME="$2"
    shift # past argument
    ;;
    -a|--alltogether)
    ALLTOGETHER="$2"
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
# 1. copy fastq files from our sever to htcf
# rsync -avFP yliang@10.20.127.3:/path/to/file .
# parallel_GNU -j 10 < fastq_rsync.sh

##############################################################
# 2. combine 2 fastq files generated in 2 different lanes (optional)
# zcat sample.Lane1.R1.fastq.gz sample.Lane2.R1.fastq.gz | gzip > sample.R1.fastq.gz
# parallel_GNU -j 10 < fastq_merge.sh
# rm *Lane*

# Example:
# rsync -avFP yliang@10.20.127.3:/taproom/data/dli/NovaSeq-S4-DT-2861687/xfer.genome.wustl.edu/gxfer1/20072858207592/CTTGGAA_S16_L003_R1_001.fastq.gz HiC_93VU147T_HPV_BRep1_TRep1_Lane1_R1.fastq.gz 
# zcat HiC_93VU147T_HPV_BRep1_TRep1_Lane1_R1.fastq.gz HiC_93VU147T_HPV_BRep1_TRep1_Lane2_R1.fastq.gz | gzip > HiC_93VU147T_HPV_BRep1_TRep1_R1.fastq.gz

##############################################################
# 3. run juicer
# (1) Get scripts from our sever
# run this manually on our server
# rsync -avFP /bar/yliang/softwares/juicer_on_htcf_myown/*sbatch yonghao@65.254.100.82:/home/yonghao/tricks/juicer_on_htcf_myown # run this manually on our server
# touch /home/yonghao/tricks/juicer_on_htcf_myown/*
# run this manually on our server

# this won't work in this direction anymore
##rsync -avFP yliang@10.20.127.3:/bar/yliang/softwares/juicer_on_htcf/ /scratch/twlab/yliang/juicer_scripts
##rsync -avFP yliang@10.20.127.3:/bar/yliang/softwares/juicer_on_htcf_myown /scratch/twlab/yliang/juicer_scripts_myown 

# (2) Run juicer on each replicate and clean up (submit job one by one based on the generation of fincin.out)
cd $INPUT_DIR
if [[ $ALLTOGETHER == "no" ]]; then
    for FILE1 in ./fastq/*${READ1EXTENSION}
    do
        xbase="$(basename $FILE1 $READ1EXTENSION)"
        echo $xbase
        mkdir $xbase $xbase/fastq
        cd $xbase/fastq
        ln -s ../../fastq/${xbase}* .
        cd ..
        echo "start $xbase juicer run"
        JUICER=run-juicer-${GENOME}_${ENZYME}.sbatch
        sbatch /home/yonghao/tricks/juicer_on_htcf_myown/$JUICER $xbase
        # sbatch /scratch/twlab/yliang/juicer_scripts/run-juicer-hg38.sbatch HiC_93VU147T_HPV_BRep1_TRep1
        while true; do
            #FILENAME=`inotifywait -m --quiet --recursive -e create,modify,close,move --format '%f' ./`
            #if [[ $FILENAME =~ fincln-[0-9]{8}.out ]]; then
            if [[ `ls ./debug/fincln-*.out > /dev/null 2> /dev/null; echo $?` -eq 0 ]]; then
                sbatch /home/yonghao/tricks/juicer_on_htcf_myown/run-cleanup.sbatch $xbase
                # sbatch /scratch/twlab/yliang/juicer_scripts/run-cleanup.sbatch HiC_93VU147T_HPV_BRep1_TRep1
                echo "$xbase juicer run finished"
                break
            fi
            #if [[ -f "./debug" ]]
        done
        cd ..
    done
else
    for FILE1 in ./fastq/*${READ1EXTENSION}
    do
        xbase="$(basename $FILE1 $READ1EXTENSION)"
        echo $xbase
        mkdir $xbase $xbase/fastq
        cd $xbase/fastq
        ln -s ../../fastq/${xbase}* .
        cd ..
        echo "start $xbase juicer run"
        JUICER=run-juicer-${GENOME}_${ENZYME}.sbatch
        sbatch /home/yonghao/tricks/juicer_on_htcf_myown/$JUICER $xbase
        cd ..
    done    
fi

# (2) Run juicer on each replicate and clean up (submit job one by one based on JOBID. However, Juicer will submit new jobs internally, which makes this not feasible)
#first_job=1
#for FILE1 in ./fastq/*${READ1EXTENSION}
#do
#    FILE2=${FILE1/R1/R2}
#    xbase="$(basename $FILE1 $READ1EXTENSION)"
#    echo $xbase
#    mkdir $xbase $xbase/fastq
#    cd $xbase/fastq
#    ln -s ../../$FILE1 .
#    ln -s ../../$FILE2 .
#    cd ..
#    if [[ $first_job == "1" ]]
#    then
#        JOBID=`sbatch --parsable /scratch/twlab/yliang/juicer_scripts/run-juicer-hg38.sbatch $xbase`
#        # sbatch --parsable /scratch/twlab/yliang/juicer_scripts/run-juicer-hg38.sbatch HiC_93VU147T_HPV_BRep1_TRep1
#        echo $JOBID
#        JOBID=`sbatch --parsable --dependency=afterok:$JOBID /scratch/twlab/yliang/juicer_scripts/run-cleanup.sbatch $xbase`
#        # sbatch --dependency=afterok:28974565 /scratch/twlab/yliang/juicer_scripts/run-cleanup.sbatch HiC_93VU147T_HPV_BRep1_TRep1
#        echo $JOBID
#        first_job=0
#    else
#        JOBID=`sbatch --parsable --dependency=afterok:$JOBID /scratch/twlab/yliang/juicer_scripts/run-juicer-hg38.sbatch $xbase`
#        echo $JOBID
#        JOBID=`sbatch --parsable --dependency=afterok:$JOBID /scratch/twlab/yliang/juicer_scripts/run-cleanup.sbatch $xbase`
#        echo $JOBID
#    fi
#    cd ..
#done

## when juicer is finished
# fincin.err is empty
# fincin.out will say you are finished
# hic.out will say everything is complete and has a time stamp at the end of the file
## after finishing generating .hic file
# clean up the folder
# sbatch /home/yonghao/tricks/juicer_on_htcf_myown/run-cleanup.batch 

# (3) Make mega map (merging replicates)
# BRep1
# --> BRep1_TRep1
#       -->aligned
#           --> merged_nodups.txt.gz
#       -->debug
# --> BRep1_TRep2

# run it under /BRep1
# sbatch /home/yonghao/tricks/juicer_on_htcf_myown/run-htcf_mega_hg38.sbatch .




