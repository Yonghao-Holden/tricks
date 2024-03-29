#!/bin/bash
# programmer: Holden Liang
# to run this script, use the wgbs_pipe_htcf_starter.sh
# this script uses bismark to process wgbs data on htcf
###################################################
# run this script on htcf to run bismark pipeline, will generate methylation level for CpG and basic qc plots including correlation heatmap across samples and methylation level distribution.
# Written by Yonghao Liang(Holden)
# for sbatch flags and how to run jobs one by one, check out websites below:
# https://slurm.schedmd.com/sbatch.html
# https://hpc.nih.gov/docs/job_dependencies.html
###################################################

#########################################################################################################################################################
## INPUT
# Arguments to bash script to decide on various parameters
INPUT_DIR="." 
READ1EXTENSION="_R1.fastq.gz"
GENOME="hg38"
SAMPLE="sample"
MAXJOBSALIGN=12 # Cores for each alignment

PBAT="no"

# default arguments
if [[ $GENOME == "hg38" ]]; then
	GENOMEDIR="/scratch/twlab/yliang/genome/hg38/bismark"
elif [[ $GENOME == "mm10" ]]; then
	GENOMEDIR="/scratch/twlab/yliang/genome/mm10/bismark"
fi

LAMBDAGENOMEDIR="/scratch/twlab/yliang/genome/lambda"
MIN_INSERT=0
MAX_INSERT=2000

fastqdir=${INPUT_DIR}/fastq
debugdir=${INPUT_DIR}/debug
trimdir=${INPUT_DIR}/trimmed
lambdadir=${INPUT_DIR}/lambda_aligned
aligndir=${INPUT_DIR}/aligned
methyldir=${INPUT_DIR}/methylcall
long_queue="general"
long_queue_time="10080"

mkdir ${debugdir} ${trimdir} ${lambdadir} ${aligndir} ${methyldir}

#########################################################################################################################################################

while [[ $# > 1 ]]
do
key="$1"
case $key in
	-i|--input)
	INPUT_DIR="$2"
	shift # past argument
	;;
	-p|--pbat)
	PBAT="$2"
	shift # past argument
	;;
	-s|--sample)
	SAMPLE="$2"
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

echo "sample name: ${SAMPLE}" >> run_info.txt
echo "fastqfile location: ${INPUT_DIR}/fastq" >> run_info.txt
echo "read 1 extension: ${READ1EXTENSION}" >> run_info.txt
echo "reference genome for bismark: ${GENOMEDIR}" >> run_info.txt

## save softwares_version information
	jid=`sbatch <<- VERSION | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/version-%j.out
	#SBATCH -e $debugdir/version-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "0_version_${SAMPLE}"

	ml fastx_toolkit bismark bwa samtools fastqc
	
	date
	echo "fastx_toolkit version:" >> softwares_version.txt
	fastx_trimmer -h >> softwares_version.txt
	echo "bismark version:" >> softwares_version.txt
	bismark -v >> softwares_version.txt
	echo "bwa version:" >> softwares_version.txt
	bwa 2>> softwares_version.txt
	echo "samtools version:" >> softwares_version.txt
	samtools --version >> softwares_version.txt
	echo "fastqc version:" >> softwares_version.txt
	fastqc -v >> softwares_version.txt
	date
VERSION`

## start processing
	jid=`sbatch <<- TRIMMING | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l 
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/trimming-%j.out
	#SBATCH -e $debugdir/trimming-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "1_trimming_${SAMPLE}"

	ml fastx_toolkit
	
	date
	find ${fastqdir} -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=\\\$(basename \\\$file) ; echo "gunzip -c "\\\$file" | fastx_trimmer -f 10 -l 150 -Q33 -z -o ${trimdir}/Trimmed_"\\\$xbase" " >> ${trimdir}/1_cutadaptcommands.txt ; done;
	find ${fastqdir} -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=\\\$(basename \\\$file) ; echo "gunzip -c "\\\${file/R1/R2}" | fastx_trimmer -f 16 -l 150 -Q33 -z -o ${trimdir}/Trimmed_"\\\${xbase/R1/R2}" " >> ${trimdir}/1_cutadaptcommands.txt ; done;
	parallel_GNU -j 2 < ${trimdir}/1_cutadaptcommands.txt
	date
TRIMMING`

dependtrimming="afterok:$jid"

	jid=`sbatch <<- FASTQC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/fastqc-%j.out
	#SBATCH -e $debugdir/fastqc-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "2_fastqc_${SAMPLE}"
	#SBATCH -d $dependtrimming    

	date
	ml fastqc
	find ${trimdir} -maxdepth 1 -name "*.gz"  | while read file ; do xbase=\\\$(basename \\\$file) ; mkdir ${trimdir}/\\\${xbase%.*}_fastqc ; echo "fastqc -o "${trimdir}/\\\${xbase%.*}"_fastqc "\\\$file"" >> ${trimdir}/1b_fastqc_commands.txt; done ;
	parallel_GNU -j 2 < ${trimdir}/1b_fastqc_commands.txt
	date
FASTQC`

	jid=`sbatch <<- ALIGNLAMBDA | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/alignlambda-%j.out
	#SBATCH -e $debugdir/alignlambda-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=128G
	#SBATCH -J "3_align_to_lambda_${SAMPLE}"
	#SBATCH -d $dependtrimming       

	date
	ml bismark bwa
	find ${fastqdir} -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=\\\$(basename \\\$file) ; echo "bismark --bowtie2 -I $MIN_INSERT -X $MAX_INSERT --multicore "$MAXJOBSALIGN" -N 1 -L 28 --score_min L,0,-0.6 -o ${lambdadir} --temp_dir ${lambdadir}/tmp --gzip --nucleotide_coverage --genome "$LAMBDAGENOMEDIR" -1 "\\\$file" -2 "\\\${file/R1/R2}" &> ${lambdadir}/"\\\$xbase"_bismark_align.log" >> ${lambdadir}/2_lambda_alignCommands.txt ; done;
	parallel_GNU -j 1 < ${lambdadir}/2_lambda_alignCommands.txt
	date
ALIGNLAMBDA`

#dependlambda="afterok:$jid"
	jid=`sbatch <<- ALIGN | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/align-%j.out
	#SBATCH -e $debugdir/align-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=128G
	#SBATCH -J "4_align_${SAMPLE}"
	#SBATCH -d $dependtrimming         

	date
	ml bismark bwa
	if [[ $PBAT == "no" ]];then
	    find ${trimdir} -name "*${READ1EXTENSION}" | while read file; do xbase=\\\$(basename \\\$file); echo "bismark --bowtie2 -I $MIN_INSERT -X $MAX_INSERT --multicore "$MAXJOBSALIGN" -N 1 -L 28 --score_min L,0,-0.6 -o ${aligndir} --temp_dir ${aligndir}/tmp --gzip --nucleotide_coverage --genome "$GENOMEDIR" -1 "\\\$file" -2 "\\\${file/R1/R2}" &> ${aligndir}/"\\\$xbase"_bismark_align.log " >> ${aligndir}/2_alignCommands.txt ; done ;
	else
	    find ${trimdir} -name "*${READ1EXTENSION}" | while read file; do xbase=\\\$(basename \\\$file); echo "bismark -q -I $MIN_INSERT -X $MAX_INSERT --parallel 2 -p 4 --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 --pbat -o ${aligndir} --nucleotide_coverage $GENOMEDIR -1 \\\$file -2 \\\${file/R1/R2} &> ${aligndir}/"\\\$xbase"_bismark_align.log " >> ${aligndir}/2_alignCommands.txt ; done ;
	fi
	parallel_GNU -j 1 < ${aligndir}/2_alignCommands.txt
	date
ALIGN`

dependalign="afterok:$jid"

	jid=`sbatch <<- DEDUP | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/dedup-%j.out
	#SBATCH -e $debugdir/dedup-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "5_dedup_${SAMPLE}"
	#SBATCH -d $dependalign      

	date
	ml bismark
	find ${aligndir} -name "*.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "deduplicate_bismark -p --output_dir ${aligndir} --bam "\\\$file" &> ${aligndir}/"\\\$xbase"_bismark_deduplication.log "  >> ${aligndir}/3_deduplicationCommands.txt ; done ;
	parallel_GNU -j 1 < ${aligndir}/3_deduplicationCommands.txt
	date
DEDUP`

dependdedup="afterok:$jid"

	jid=`sbatch <<- methylcall | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l 
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/methylextract-%j.out
	#SBATCH -e $debugdir/methylextract-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=128G
	#SBATCH -J "6_methylcall_${SAMPLE}"
	#SBATCH -d $dependdedup 

	date
	ml bismark

	find ${aligndir} -name "*.deduplicated.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report --bedgraph --cytosine_report -o ${methyldir} --gzip --multicore "$MAXJOBSALIGN" "\\\$file" &> ${methyldir}/"\\\$xbase"_bismark_methylcall.log --genome_folder $GENOMEDIR" >> ${methyldir}/4_extractorCommands.txt ; done ;
	parallel_GNU -j 1 < ${methyldir}/4_extractorCommands.txt
	date
methylcall`

dependmethyl="afterok:$jid"

	jid=`sbatch <<- final | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/final-%j.out
	#SBATCH -e $debugdir/final-%j.err 
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "7_final_${SAMPLE}"
	#SBATCH -d $dependmethyl  

	date
	echo "yes, it's done!" > $debugdir/done-${SAMPLE}.out
	date
final`
