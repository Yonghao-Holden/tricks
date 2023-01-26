#!/bin/bash
# programmer: Holden (Yonghao)
# This is a simple wrapper script that will allow for running the SQuIRE pipeline

# Arguments to bash script to decide on various parameters
MAXJOBS=6 #Cores per CPU
STARJOB=2 #Cores per alignment
FASTQFOLDER="." # folder that contains all fastq files # use absolute path
READ1EXTENSION="R1.fastq.gz"
FETCHFOLDER="/bar/yliang/softwares/SQuIRE/squire_fetch"
CLEANFOLDER="/bar/yliang/softwares/SQuIRE/squire_clean"
GTF="/bar/yliang/genomes/private/gencode.v29.annotation.gtf"
BUILD="hg38"
READLENGTH=75
STRNAD=2 # juheon tested this, need to check in stringtie mannual to see which one is the correct one.

while [[ $# > 1 ]]
do
key="$1"
case $key in
    -j|--maxjobs)
    MAXJOBS="$2"
    shift # past argument
    ;;
    --starjob)
    STARJOB="$2"
    shift # past argument
    ;;
    --fastqfolder)
    FASTQFOLDER="$2"
    shift # past argument
    ;;
	-r|--read1extension)
    READ1EXTENSION="$2"
    shift # past argument
    ;;
    --fetchfolder)
    FETCHFOLDER="$2"
    shift # past argument
    ;;
    --cleanfolder)
    CLEANFOLDER="$2"
    shift # past argument
    ;;
    --gtf)
    GTF="$2"
    shift # past argument
    ;;
    --build)
    BUILD="$2"
    shift # past argument
    ;;
    --readlength)
    READLENGTH="$2"
    shift # past argument
    ;;
    --strand)
    STRAND="$2"
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

# conda create --name squire --override-channels -c iuc -c bioconda -c conda-forge -c defaults -c r python=2.7.13 bioconductor-deseq2=1.16.1 
# r-base=3.4.1 r-pheatmap bioconductor-vsn bioconductor-biocparallel=1.12.0 r-ggrepel star=2.5.3a bedtools=2.25.0 samtools=1.1 stringtie=1.3.3 igvtools=2.3.93 
# ucsc-genepredtobed ucsc-gtftogenepred ucsc-genepredtogtf ucsc-bedgraphtobigwig r-hexbin

echo "parallel job limit: ${MAXJOBS}"

ml cutadapt FastQC

echo "STAR version:" >> softwares_version.txt
STAR --version >> softwares_version.txt
echo "bedtools version:" >> softwares_version.txt
bedtools --version >> softwares_version.txt
echo "samtools version:" >> softwares_version.txt
samtools --version >> softwares_version.txt
echo "R version:" >> softwares_version.txt
R --version >> softwares_version.txt
echo "Version information for other softwres can be seen in the conda environment setting" >> softwares_version.txt

##########################################################################################
mkdir trimmed
cd trimmed

find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o Trimed_"$xbase" -p Trimed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 
echo "(1/4) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt

##########################################################################################
find ../trimmed -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc "$file >> 2_fastqc_commands.txt; done ;
echo "(2/4) FastQC run"
parallel_GNU -j $MAXJOBS < 2_fastqc_commands.txt

##########################################################################################
cd ..
source activate squire

# Squire set up
# only needs to be run once
# /bar/yliang/softwares/SQuIRE
#squire Fetch -b hg38 -f -c -r -g -x -p 10 
#squire Clean -b hg38

# Mapping
find ./trimmed -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file); echo "squire Map -1 $file -2 ${file/R1/R2} -o squire_map -f $FETCHFOLDER -r $READLENGTH -n $xbase -b $BUILD -g $GTF -p $STARJOB" >> 3_squire_map.txt ; done ;
# find ./trimmed -name "*R1.fastq.gz" | while read file ; do xbase=$(basename $file); echo "squire Map -1 $file -2 ${file/R1/R2} -o squire_map -f /bar/yliang/softwares/SQuIRE/squire_fetch -r 75 -n $xbase -b hg38 -g /bar/yliang/genomes/private/gencode.v29.annotation.gtf -p 2" >> 3_squire_map.txt ; done ;
echo "(3/4) Mapping Reads"
parallel_GNU -j $MAXJOBS < 3_squire_map.txt

# TE Counting
find ./trimmed -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file); echo "squire Count -m squire_map -c $CLEANFOLDER -o squire_count -f $FETCHFOLDER -r $READLENGTH -n $xbase -b $BUILD -p 10 -s $STRAND" >> 4_squire_count.txt; done; 
# find ./trimmed -name "*R1.fastq.gz" | while read file ; do xbase=$(basename $file); echo "squire Count -m squire_map -c /bar/yliang/softwares/SQuIRE/squire_clean -o squire_count -f /bar/yliang/softwares/SQuIRE/squire_fetch -r 75 -n $xbase -b hg38 -p 10 -s 2" >> 4_squire_count.txt; done; 
echo "(4/4) Counting TEs"
parallel_GNU -j $MAXJOBS < 4_squire_count.txt 

# DE-Seq2
# squire Call -1 Trimed_RNA_kera_3O2_BRep1.R1.fastq.gz -2 Trimed_JHU006_NEG_R1.fastq.gz -A kera -B NEG -o kera_vs_JHU006_NEG_call -p 10 -f pdf -i squire_count 

conda deactivate

##########################################################################################
# plot good looking volcano plot
# Rscript /bar/yliang/tricks/squire_2.R -wd $PWD -data all -output AvsB_volcano_TE.pdf

