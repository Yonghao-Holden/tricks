#!/bin/bash -e

# programmer: Holden
# This is a pipline to do the initial steps with wgbs data up to assembly
########################################################################
# Arguments to bash script to decide on various parameters
MAXJOBS=6 #Cores per CPU
MAXJOBSALIGN=4 # Cores for each alignment
FASTQFOLDER="./fastq"
READ1EXTENSION="R1.fastq.gz"
########################################################################
GENOMEDIR="/bar/yliang/genomes/public/hg38/bismark"
LAMBDAGENOMEDIR="/bar/yliang/genomes/private/lambda_bismark"
MIN_INSERT=0
MAX_INSERT=2000

while [[ $# > 1 ]]
do
key="$1"
case $key in
    -j|--maxjobs)
    MAXJOBS="$2"
    shift # past argument
    ;;
    -aj|--maxjobsalign)
    MAXJOBSALIGN="$2"
    shift # past argument
    ;;
    -f|--fastqfolder)
    FASTQFOLDER="$2"
    shift # past argument
    ;;
    -r|--read1extension)
    READ1EXTENSION="$2"
    shift # past argument
    ;;
    -g|--genomedir)
    GENOMEDIR="$2"
    shift # past argument
    ;;
    -lg|--lambdagenomedir)
    LAMBDAGENOMEDIR="$2"
    shift # past argument
    ;;
    -mini|--mininsert)
    MIN_INSERT="$2"
    shift # past argument
    ;;
    -maxi|--maxinsert)
    MAX_INSERT="$2"
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

########################################################################
echo "parallel job limit: ${MAXJOBS}"
echo "parallel job limit for alignment: ${MAXJOBSALIGN}"
echo "fastqfile location: ${FASTQFOLDER}"
echo "read 1 extension: ${READ1EXTENSION}"
echo "reference genome for bismark: ${GENOMEDIR}"

ml bismark bwa FastQC

## save softwares_version information
cd $FASTQFOLDER/..
echo "fastx_toolkit version:" >> softwares_version.txt
fastx_trimmer -h >> softwares_version.txt
echo "bismark version:" >> softwares_version.txt
bismark -v >> softwares_version.txt
echo "bwa version:" >> softwares_version.txt
bwa 2>> softwares_version.txt
echo "samtools version:" >> softwares_version.txt
samtools 2>> softwares_version.txt
echo "fastqc version:" >> softwares_version.txt
fastqc -v >> softwares_version.txt
########################################################################
## 1. cut adaptor
mkdir trimmed
cd trimmed
#find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o Trimed_"$xbase" -p Trimed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 
find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "gunzip -c "$file" | fastx_trimmer -f 10 -l 150 -Q33 -z -o Trimmed_"$xbase" " >> 1_cutadaptcommands.txt ; done;
find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "gunzip -c "${file/R1/R2}" | fastx_trimmer -f 16 -l 150 -Q33 -z -o Trimmed_"${xbase/R1/R2}" " >> 1_cutadaptcommands.txt ; done;
echo "(1/5) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt

## 1b. fastqc
find ./ -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc "$file"" >> 1b_fastqc_commands.txt; done ;
echo "(1b/5) FastQC run"
parallel_GNU -j $MAXJOBS < 1b_fastqc_commands.txt

########################################################################
## 2a. align to lambda genome
cd ..
mkdir lambda_aligned
cd lambda_aligned
#bismark --bowtie2 --multicore 8 -N 1 --score_min L,0,-0.2 -o ./lambda_aligned /scratch/twlab/hlee/genomes/lambda -1 ./trimmedreads/fastx_trim/${part1}'_R1_fastxtrimmed.fq.gz' -2 ./trimmedreads/fastx_trim/${part1}'_R2_fastxtrimmed.fq.gz'
find ../trimmed -maxdepth 1 -name "Trimmed_*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "bismark --bowtie2 --multicore "$MAXJOBSALIGN" -N 1 --score_min L,0,-0.2 -o ./ --genome "$LAMBDAGENOMEDIR" -1 "$file" -2 "${file/R1/R2}" &> ./"$xbase"_bismark_align.log" >> 2a_lambda_alignCommands.txt ; done;
echo "(2a/5) Aligning Reads to lambda genome"
parallel_GNU -j $MAXJOBS < 2a_lambda_alignCommands.txt

########################################################################
## 2b. align to human genome
cd ..
mkdir aligned
cd aligned
#bismark -I $MIN_INSERT -X $MAX_INSERT --parallel $CPU -o $outdir --temp_dir ./tmp --gzip --nucleotide_coverage $genome_dir -1 $inputR1 -2 $inputR2 &>/scratch/twlab/jjang/TE_antigen/WGBS/aligned/${part1}'_bismark.log'
find ../trimmed -name "Trimmed_*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "bismark --bowtie2 -I $MIN_INSERT -X $MAX_INSERT --multicore "$MAXJOBSALIGN" -o ./ --temp_dir ./"$xbase"_tmp --gzip --nucleotide_coverage --genome "$GENOMEDIR" -1 "$file" -2 "${file/R1/R2}" &> ./"$xbase"_bismark_align.log " >> 2b_alignCommands.txt ; done ;
echo "(2b/5) Aligning Reads to human genome"
parallel_GNU -j $MAXJOBS < 2b_alignCommands.txt

########################################################################
## 3. call methylation
#bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report -o $outdir2 --gzip --parallel $CPU $bam_dedup &>$log_methx
find ./ -name "*.bam" | while read file; do xbase=$(basename $file); echo "bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report --gzip --bedGraph --multicore "$MAXJOBSALIGN" "$file" &> "$xbase"_bismark_methylcall.log " >> 3_extractorCommands.txt ; done ;
echo "(3/5) Methylation Extraction"
parallel_GNU -j $MAXJOBS < 3_extractorCommands.txt

########################################################################
## 4. bismark report and summary and multiQC
bismark2report
bismark2summary

ml samtools
find ./ -name "*pe.bam" | while read file; do xbase=${file/.bam/}; echo "samtools sort $file -o $xbase.sorted.bam; samtools index $xbase.sorted.bam" >> 4_sort_index_bam_Commnads.txt; done;
parallel_GNU -j $MAXJOBS < 4_sort_index_bam_Commnads.txt

cd ..
multiqc .
