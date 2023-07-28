#!/bin/bash
# programmer: Nakul
# This is a pipline to do the initial steps with RNA-seq data up to assembly and then the main pipeline that we have can be run

# Arguments to bash script to decide on various parameters
MAXJOBS=6 #Cores per CPU
MAXJOBSSTAR=2 #STAR requires a lot of resources so it is not revommended to have more than 2 alignment running at the same time. 
FASTQFOLDER="./"
READ1EXTENSION="R1.fastq.gz"

while [[ $# > 1 ]]
do
key="$1"
case $key in
    -j|--maxjobs)
    MAXJOBS="$2"
    shift # past argument
    ;;
    -s|--maxstarjobs)
    MAXJOBSSTAR="$2"
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
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

echo "parallel job limit: ${MAXJOBS}"
echo "parallel job limit STAR: ${MAXJOBSSTAR}"
echo "fastqfile location: ${FASTQFOLDER}"
echo "read 1 extension: ${READ1EXTENSION}"

ml cutadapt STAR FastQC stringtie python2 multiqc picard

mkdir trimmed
cd trimmed

find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o Trimed_"$xbase" -p Trimed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 

echo "(1/5) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt

find ../trimmed -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc "$file >> fastqc_commands.txt; done ;

echo "(1b/5) FastQC run"
parallel_GNU -j $MAXJOBS < fastqc_commands.txt

cd ..
mkdir aligned
cd aligned

find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "STAR  --runMode alignReads  --runThreadN 4  --genomeDir ~/reference/STAR_index_hg38_gencodeV22/ --readFilesIn "$file" "${file/R1/R2}" --readFilesCommand zcat --outFileNamePrefix "${xbase%.*}" --outSAMtype BAM   SortedByCoordinate   --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMheaderHD @HD VN:1.4 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33  --alignIntronMax 500000  --alignMatesGapMax 1000000 --twopassMode Basic" >> 2_alignCommands.txt ; done ;

echo "(2/5) Aligning Reads with STAR"
parallel_GNU -j $MAXJOBSSTAR < 2_alignCommands.txt

ls *bam | while read file ; do samtools index $file ; samtools idxstats $file > ${file}_idxstats ; done ;

cd ..
mkdir rseqc
cd rseqc

find ../aligned -name "*bam" | while read file ; do xbase=$(basename $file) ; echo "geneBody_coverage.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38.HouseKeepingGenes.bed" >> 3_rseqQCcommands.txt ; echo "read_distribution.py -i "$file" -r ~/reference/rseqc/hg38_Gencode_V28.bed > "${xbase%.*}".readdistribution.txt" >> 3_rseqQCcommands.txt ; echo "junction_saturation.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38_Gencode_V28.bed" >> 3_rseqQCcommands.txt ; echo "bam2wig.py -s "~/reference/mapability/hg38.chrom.sizes" -i "$file" -o "${xbase%.*}" -u" >> 3_rseqQCcommands.txt; done ;

echo "(3/5) Quality Control and Visualization"
parallel_GNU -j $MAXJOBS < 3_rseqQCcommands.txt

find ../aligned -name "*bam" | while read file ; do echo "java -jar $PICARD EstimateLibraryComplexity I="$file" O="${file%.*}"_duplication_stats.txt" >> estimateComplexity.txt ; done ;

echo "(3.5/5) Estimate Complexity"
parallel_GNU -j $MAXJOBS < estimateComplexity.txt

cd ..
mkdir assembled
cd assembled

find ../aligned -name "*bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -q 255 -h "$file" | stringtie - -o "${xbase%.*}".gtf -p 4 -m 100 -c 1 --rf" >> 4_assembleCommands.txt ; done ;

echo "(4/5) Quality Control and Visualization"
parallel_GNU -j $MAXJOBS < 4_assembleCommands.txt

cd ..
echo "(5/5) multiQC"

multiqc .

