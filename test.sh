#!/bin/bash
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, none). Optional arguments are the queue for the 
# alignment, description for stats file, 
# stage to relaunch at, paths to various files if needed,
# chunk size, path to scripts directory, and the top-level directory (default 
# current directory). In lieu of setting the genome ID, you can instead set the
# reference sequence and the chrom.sizes file path, but the directory 
# containing the reference sequence must also contain the BWA index files.
#
# Splits the fastq files, creates jobs to align them, creates merge jobs that
# wait for the alignment to finish, and creates a final merge job.
#
# Also creates "cleanup" jobs that at each stage, deletes jobs off the cluster
# if any one of them fails.
#
# If all is successful, takes the final merged file, removes name duplicates,
# removes PCR duplicates, and creates the hic job and stats job.  Final
# product will be hic file and stats file in the aligned directory.
#                                                                       
# [topDir]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "R" in the appropriate files, i.e. *R*.fastq
# From the top-level directory, the following two directories are created:
#                                                                              
# [topDir]/splits  - Where to write the scratch split files (fastq files and
#                    intermediate SAM files). This can be deleted after 
#                    execution.
# [topDir]/aligned - Where to write the final output files.
#
# The following globals should be set correctly before proceeding:
#
# splitsize - The number of lines that each split fastq should contain. Larger
#             means fewer files and longer overall, but too small means there
#             are so many jobs that the cluster won't run them. This can be
#             set with the -C command as well
# read1str  - portion of fastq filename that indicates this is the "read 1"
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly,
#             script will not work. The error will often manifest itself
#             through a "*" in the name because the wildcard was not able to
#             match any files with the read1str.   
# Juicer version 1.6
shopt -s extglob
juicer_version="1.6"
## Set the following variables to work with your system

# BCK: added 2020-08-20 so slurm doesn't complain
unset SLURM_MEM_PER_CPU

# Aiden Lab specific check
isRice=$(host $(hostname) | awk '{if ($1~/rice/) {print 1} else {print 0}; exit}') #'
isBCM=$(host $(hostname) | awk '{if ($1~/bcm/) {print 1} else {print 0}; exit}') #'
isVoltron=0
## path additionals, make sure paths are correct for your system
## use cluster load commands
if [ $isRice -eq 1 ] 
then
	export PATH=$HOME/bin:$PATH
	isNots=$(host $(hostname) | awk '{if ($1~/nots/){print 1}else {print 0}}') #'
	if [ $isNots -eq 1 ]
	then
	load_bwa="module load  GCCcore/7.3.0 BWA/0.7.17"
	load_java="module load Java/1.8.0_162" 
	load_gpu="module load gcccuda/2016a;module load CUDA/8.0.44;" 
	else
	load_bwa="export PATH=/home/ncd4/bwa:$PATH"
	load_java="module load Java/8.0.3.22" 
	load_gpu="module load gcccuda/2016a;module load CUDA/8.0.54;" 
	fi
	
	# Juicer directory, contains scripts/, references/, and restriction_sites/
	# can also be set in options via -D
	juiceDir="/projects/ea14/juicer" ### RICE
	# default queue, can also be set in options via -q
	queue="commons"
	queue_time="24:00:00"
	# default long queue, can also be set in options via -l
	long_queue="commons"
	long_queue_time="24:00:00"
elif [ $isBCM -eq 1 ]
then    
	# Juicer directory, contains scripts/, references/, and restriction_sites/
	# can also be set in options via -D
	juiceDir="/storage/aiden/juicer/"
	# default queue, can also be set in options via -q
	queue="mhgcp"
	queue_time="1200"
	# default long queue, can also be set in options via -l
	long_queue="mhgcp"
	long_queue_time="3600"
else
	isVoltron=1
	#export PATH=/gpfs0/biobuild/biobuilds-2016.11/bin:$PATH 
	# unset MALLOC_ARENA_MAX # on IBM platform this parameter has significant speed efect but may result in memory leaks
	unset MALLOC_ARENA_MAX
	load_bwa="module load bwa" 
	load_java="module load java" 
	load_gpu="module load cuda/8.0.61; CUDA_VISIBLE_DEVICES=0,1,2,3"
	#call_gem="/gpfs0/work/neva/gem3-mapper/bin/gem-mapper --3c"
	# Juicer directory, contains scripts/, references/, and restriction_sites/
	# can also be set in options via -D
	juiceDir="/home/yonghao/softwares/juicer/"
	# default queue, can also be set in options
	queue="general"
	queue_time="1440"
	# default long queue, can also be set in options
	long_queue="general"
	long_queue_time="10080"
fi

# size to split fastqs. adjust to match your needs. 4000000=1M reads per split
# can also be changed via the -C flag
splitsize=90000000

# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 

# unique name for jobs in this run
groupname="a$(date +%s)"

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
# default: not set
site="MboI"
# genome ID, default to human, can also be set in options
genomeID="hg38"
# description, default empty
about=""
# do not include fragment delimited maps by default
nofrag=0
# use wobble for dedupping by default (not just exact matches)
justexact=0

## Read arguments                                                     
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-C chunk size] [-D Juicer scripts directory]\n                 [-Q queue time limit] [-L long queue time limit] [-b ligation] [-t threads]\n                 [-A account name] [-e] [-h] [-f] [-j]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \n  \"$genomeID\"); alternatively, it can be defined using the -z command"
dirHelp="* [topDir] is the top level directory (default\n  \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
queueHelp="* [queue] is the queue for running alignments (default \"$queue\")"
longQueueHelp="* [long queue] is the queue for running longer jobs such as the hic file\n  creation (default \"$long_queue\")"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" \n  (default \"$site\")"
aboutHelp="* [about]: enter description of experiment, enclosed in single quotes"
stageHelp="* [stage]: must be one of \"chimeric\", \"merge\", \"dedup\", \"final\", \"postproc\", or \"early\".\n    -Use \"merge\" when alignment has finished but the merged_sort file has not\n     yet been created.\n    -Use \"dedup\" when the files have been merged into merged_sort but\n     merged_nodups has not yet been created.\n    -Use \"final\" when the reads have been deduped into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the hic files\n    Can also use -e flag to exit early"
pathHelp="* [chrom.sizes path]: enter path for chrom.sizes file"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
chunkHelp="* [chunk size]: number of lines in split files, must be multiple of 4\n  (default ${splitsize}, which equals $(awk -v ss=${splitsize} 'BEGIN{print ss/4000000}') million reads)"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
refSeqHelp="* [reference genome file]: enter path for reference sequence file, BWA index\n  files must be in same directory"
queueTimeHelp="* [queue time limit]: time limit for queue, i.e. -W 12:00 is 12 hours\n  (default ${queue_time})"
longQueueTimeHelp="* [long queue time limit]: time limit for long queue, i.e. -W 168:00 is one week\n  (default ${long_queue_time})"
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
threadsHelp="* [threads]: number of threads when running BWA alignment"
userHelp="* [account name]: user account name on cluster"
excludeHelp="* -f: include fragment-delimited maps in hic file creation"
justHelp="* -j: just exact duplicates excluded at dedupping step"
earlyexitHelp="* -e: Use for an early exit, before the final creation of the hic files"
gemHelp="* -c: use GEM3 as aligner"
helpHelp="* -h: print this help and exit"


printHelpAndExit() {
	echo -e "$usageHelp"
	echo -e "$genomeHelp"
	echo -e "$dirHelp"
	echo -e "$queueHelp"
	echo -e "$longQueueHelp"
	echo -e "$siteHelp"
	echo -e "$aboutHelp"
	echo -e "$stageHelp"
	echo -e "$pathHelp"
	echo -e "$siteFileHelp"
	echo -e "$refSeqHelp"
	echo -e "$chunkHelp"
	echo -e "$scriptDirHelp"
	echo -e "$queueTimeHelp"
	echo -e "$longQueueTimeHelp"
	echo -e "$ligationHelp"
	echo -e "$threadsHelp"
	echo -e "$userHelp"
	echo -e "$earlyexitHelp"
	echo -e "$gemHelp"
	echo -e "$justHelp"
	echo "$excludeHelp"
	echo "$helpHelp"
	exit "$1"
}

while getopts "d:g:a:hq:s:p:l:y:z:S:C:D:Q:L:b:A:t:jfec" opt; do
	case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	l) long_queue=$OPTARG ;;
	q) queue=$OPTARG ;;
	s) site=$OPTARG ;;
	a) about=$OPTARG ;;
	p) genomePath=$OPTARG ;;  
	y) site_file=$OPTARG ;;
	z) refSeq=$OPTARG ;;
	S) stage=$OPTARG ;;
	C) splitsize=$OPTARG; splitme=1 ;;
	D) juiceDir=$OPTARG ;;
	Q) queue_time=$OPTARG ;;
	L) long_queue_time=$OPTARG ;;
	f) nofrag=0 ;;
	b) ligation=$OPTARG ;;
	t) threads=$OPTARG ;;
	A) user=$OPTARG ;;
	j) justexact=1 ;;
	e) earlyexit=1 ;;
	c) gemmapper=1 ;;
	[?]) printHelpAndExit 1;;
	esac
done

if [ ! -z "$stage" ]
then
	case $stage in
	chimeric) chimeric=1 ;;
		merge) merge=1 ;;
		dedup) dedup=1 ;;
		early) earlyexit=1 ;;
		final) final=1 ;;
	postproc) postproc=1 ;; 
		*)  echo "$usageHelp"
		echo "$stageHelp"
		exit 1
	esac
fi

## Set reference sequence based on genome ID
if [ -z "$refSeq" ]
then 
	case $genomeID in
	mm9)	refSeq="${juiceDir}/references/Mus_musculus_assembly9_norandom.fasta";;
	mm10)	refSeq="${juiceDir}/references/Mus_musculus_assembly10/v0/Mus_musculus_assembly10.fasta";;
	hg38)	refSeq="${juiceDir}/references/hg38/hg38.fa";;
	hg19)	refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
	hg18)	refSeq="${juiceDir}/references/hg18.fasta";;
	*)	echo "$usageHelp"
		echo "$genomeHelp"
		exit 1
	esac
else
	## Reference sequence passed in, so genomePath must be set for the .hic 
	## file to be properly created
	if [[ -z "$genomePath" ]] && [[ -z $earlyexit ]]
	then
		echo "***! You must define a chrom.sizes file or a standard genome ID via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq; you may use \"-p hg19\" or other standard genomes";
		exit 1;
	fi
fi

## Check that refSeq exists 
if [ ! -e "$refSeq" ]; then
	echo "***! Reference sequence $refSeq does not exist";
	exit 1;
fi

## Check that index for refSeq exists
if [[ ! -e "${refSeq}.bwt" ]] && [[ -z $gemmapper ]]
then
	echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run bwa index on this file before running juicer.";
	exit 1;
elif [[ -n $gemmapper ]] && [[ ! -e "${refSeq%.*}.gem" ]]
then
	echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run gem index on this file before running juicer.";
	exit 1;
fi

## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]; then
	case $site in
	HindIII) ligation="AAGCTAGCTT";;
	MseI)  ligation="TTATAA";;
	DpnII) ligation="GATCGATC";;
	MboI) ligation="GATCGATC";;
		NcoI) ligation="CCATGCATGG";;
	Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
	none) ligation="XXXX";;
	*)  ligation="XXXX"
		echo "$site not listed as recognized enzyme."
		echo "Ligation junction is undefined"
	esac
fi

if [[ -n $gemmapper ]] 
then
	if [[ "$site" == "none" ]]
	then
	re=""
	else
	re=$(echo $site | tr "+" "\n" | awk '{str=str" --restriction-enzyme "$1}END{print str}')
	fi
fi

## If DNAse-type experiment, no fragment maps; or way to get around site file
if [[ "$site" == "none" ]] 
then
	nofrag=1;
fi

if [ -z "$site_file" ]
then
	site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [[ ! -e "$site_file" ]] && [[ "$site" != "none" ]] &&  [[ ! "$site_file" =~ "none" ]]
then
	echo "***! $site_file does not exist. It must be created before running this script."
	exit 1
elif [[ "$site" != "none" ]] && [[ ! "$site_file" =~ "none" ]]
then
	echo  "Using $site_file as site file"
fi

## Set threads for sending appropriate parameters to cluster and string for BWA call
if [ -z "$threads" ]
then
	# default is 8 threads; may need to adjust
	if [ $isRice -eq 1 ]
	then
	threads=8
	threadstring="-t $threads"
	elif [ $isBCM -eq 1 ]
	then
	threads=1
	threadstring="-t $threads"
	else
	threads=16 # VOLTRON; may need to make this separate and specific for Voltron
	## On voltron with 8 thread per core Power8 CPU bwa can use more threads
	threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"
	fi
else
	if [ $isVoltron -eq 1 ]
	then
	threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"
	else
	threadstring="-t $threads"
	fi
fi

alloc_mem=$(($threads * 5000))

if [ $alloc_mem -gt 40000 ]
then
	alloc_mem=40000
fi

if [ $isBCM -eq 1 ] || [ $isRice -eq 1 ]
then
	alloc_mem=50000
fi

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"
debugdir=${topDir}"/debug"

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$topDir/fastq" ]; then
	echo "Directory \"$topDir/fastq\" does not exist."
	echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
	echo "Type \"juicer.sh -h\" for help"
	exit 1
else 
	if stat -t ${fastqdir} >/dev/null 2>&1
	then
	echo "(-: Looking for fastq files...fastq files exist"
	else
	if [ ! -d "$splitdir" ]; then 
		echo "***! Failed to find any files matching ${fastqdir}"
		echo "***! Type \"juicer.sh -h \" for help"
		exit 1		
	fi
	fi
fi

## Create output directory, only if not in postproc, dedup or final stages
if [[ -d "$outputdir" && -z "$final" && -z "$dedup" && -z "$postproc" ]] 
then
	echo "***! Move or remove directory \"$outputdir\" before proceeding."
	echo "***! Type \"juicer.sh -h \" for help"
	exit 1			
else
	if [[ -z "$final" && -z "$dedup" && -z "$postproc" ]]; then
		mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
	fi
fi

## Create split directory
if [ -d "$splitdir" ]; then
	splitdirexists=1
else
	mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ]; then
	mkdir "$tmpdir"
	chmod 777 "$tmpdir"
fi

## Create output directory, used for reporting commands output
if [ ! -d "$debugdir" ]; then
	mkdir "$debugdir"
	chmod 777 "$debugdir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline
# If chunk size sent in, split. Otherwise check size before splitting
if [ $isVoltron -ne 1 ]
then
	if [ -z $splitme ]
	then
	fastqsize=$(ls -lL  ${fastqdir} | awk '{sum+=$5}END{print sum}')
	if [ "$fastqsize" -gt "2592410750" ]
	then
		splitme=1
	fi
	fi
fi

testname=$(ls -l ${fastqdir} | awk 'NR==1{print $9}')
if [ "${testname: -3}" == ".gz" ]
then
	read1=${splitdir}"/*${read1str}*.fastq.gz"
	gzipped=1
else
	read1=${splitdir}"/*${read1str}*.fastq"
fi

if [ -z "$user" ]
then
	userstring=""
else
	userstring="#SBATCH -A $user"
fi

# Add header containing command executed and timestamp:
jid=`sbatch <<- HEADER | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l 
		$userstring
	#SBATCH -p $queue
	#SBATCH -t 2
	#SBATCH -c 1
	#SBATCH -o $debugdir/head-%j.out
	#SBATCH -e $debugdir/head-%j.err
	#SBATCH -J "${groupname}_cmd"
	date
	${load_bwa}
	${load_java}
	${load_awk}

	echo "header"

HEADER`

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
    samtools 2>> softwares_version.txt
    echo "fastqc version:" >> softwares_version.txt
    fastqc -v >> softwares_version.txt
    date
VERSION`

