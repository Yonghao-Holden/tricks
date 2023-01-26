#!/bin/bash
#!/usr/bin/awk
# programmer: Nakul (adapted by Holden)
# Run it under the folder `<ATAC folder>/`
#####################################################################
# This is a pipline to perform initial analysis on RNA-seq data including de novo assemble and gene expression quantification, which will be ready to use for TEProf2 and DESeq2 directly, respectively.
#####################################################################
# Update log:
# 1. version1 used inconsistent version of gencode for alignment and gene quantification. fix this in version2
#####################################################################
# default settings
MAXJOBS=6 #Cores per CPU
MAXJOBSSTAR=2 #STAR requires a lot of resources so it is not revommended to have more than 2 alignment running at the same time. 
FASTQFOLDER="./"
READ1EXTENSION="R1.fastq.gz"
QUANTIFICATION="yes"

# Arguments to bash script to decide on various parameters
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
    -q|--quantification)
    QUATIFICATION="$2"
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

cd $FASTQFOLDER/..

echo "parallel job limit: ${MAXJOBS}"
echo "parallel job limit STAR: ${MAXJOBSSTAR}"
echo "fastqfile location: ${FASTQFOLDER}"
echo "read 1 extension: ${READ1EXTENSION}"

ml cutadapt STAR FastQC stringtie python2 picard samtools

## save softwares_version information
echo "cutadapt version:" >> softwares_version.txt
cutadapt --version >> softwares_version.txt
echo "STAR version:" >> softwares_version.txt
STAR --version >> softwares_version.txt
echo "FastQC version:" >> softwares_version.txt
fastqc --version >> softwares_version.txt
echo "stringtie version:" >> softwares_version.txt
stringtie --version >> softwares_version.txt
echo "multiQC version:" >> softwares_version.txt
multiqc --version >> softwares_version.txt
echo "picard_EstimateLibraryComplexity version:" >> softwares_version.txt
java -jar $PICARD EstimateLibraryComplexity --version 2>> softwares_version.txt
echo "samtools version:" >> softwares_version.txt
samtools --version >> softwares_version.txt
echo "featureCounts version:" >> softwares_version.txt
featureCounts -v 2>> softwares_version.txt

mkdir trimmed
cd trimmed

find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o Trimmed_"$xbase" -p Trimmed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 
echo "(1/10) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt

find ../trimmed -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc "$file >> 2_fastqc_commands.txt; done ;
echo "(2/10) FastQC run"
parallel_GNU -j $MAXJOBS < 2_fastqc_commands.txt

cd ..
mkdir aligned
cd aligned

find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "STAR  --runMode alignReads  --runThreadN 4  --genomeDir /bar/yliang/genomes/private/STAR_index_hg38_gencodeV36/ --readFilesIn "$file" "${file/R1/R2}" --readFilesCommand zcat --outFileNamePrefix "${xbase%.*}" --outSAMtype BAM   SortedByCoordinate   --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMheaderHD @HD VN:1.4 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33  --alignIntronMax 500000  --alignMatesGapMax 1000000 --twopassMode Basic" >> 3_alignCommands.txt ; done ;
echo "(3/10) Aligning Reads with STAR"
parallel_GNU -j $MAXJOBSSTAR < 3_alignCommands.txt

find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "samtools index $file ; samtools idxstats $file > ${file}_idxstats ; samtools sort $file -O bam -@ 3 -n -T ${file}_temp > ${file/sortedByCoord/sortedByReadname} " >> 4_sortingCommands.txt; done;
echo "(4/10) Sort bam files by readname"
parallel_GNU -j $MAXJOBS < 4_sortingCommands.txt

find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "/bar/yliang/anaconda3/bin/bamCoverage -b $file -o ${file/bam/reverse.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand forward --scaleFactor -1 2> 5_bam2bigwigCommand.log" >> 5_bam2bigwigCommand.txt; done;
find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "/bar/yliang/anaconda3/bin/bamCoverage -b $file -o ${file/bam/forward.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand reverse 2> 5_bam2bigwigCommand.log" >> 5_bam2bigwigCommand.txt; done;
echo "(5/10) Convert bam to bigwig"
parallel_GNU -j $MAXJOBS < 5_bam2bigwigCommand.txt

cd ..
mkdir rseqc
cd rseqc

find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "geneBody_coverage.py -i "$file" -o "${xbase%.*}" -r /bar/yliang/genomes/private/rseqc/hg38.HouseKeepingGenes.bed" >> 6_rseqQCcommands.txt ; echo "read_distribution.py -i "$file" -r /bar/yliang/genomes/private/rseqc/hg38_Gencode_V28.bed > "${xbase%.*}".readdistribution.txt" >> 6_rseqQCcommands.txt ; echo "junction_saturation.py -i "$file" -o "${xbase%.*}" -r /bar/yliang/genomes/private/rseqc/hg38_Gencode_V28.bed" >> 6_rseqQCcommands.txt ; echo "bam2wig.py -s "/bar/yliang/genomes/private/hg38_all_chromosomes.size" -i "$file" -o "${xbase%.*}" -u" >> 6_rseqQCcommands.txt; done ;
echo "(6/10) Quality Control and Visualization"
parallel_GNU -j $MAXJOBS < 6_rseqQCcommands.txt

find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do echo "java -jar \$PICARD EstimateLibraryComplexity I="$file" O="${file%.*}"_duplication_stats.txt" >> 7_estimateComplexity.txt ; done ;
echo "(7/10) Estimate Complexity"
parallel_GNU -j $MAXJOBSSTAR < 7_estimateComplexity.txt

cd ..
mkdir assembled
cd assembled

find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -q 255 -h "$file" | stringtie - -o "${xbase%.*}".gtf -p 4 -m 100 -c 1 --fr" >> 8_assembleCommands.txt ; done ;
echo "(8/10) Quality Control and Visualization"
parallel_GNU -j $MAXJOBS < 8_assembleCommands.txt

cd ..
echo "(9/10) multiQC"

/bar/yliang/anaconda3/bin/multiqc .

mv multiqc_report.html $(basename "$PWD").multiqc_report.html

if [[ $QUANTIFICATION == "yes" ]]
then
    mkdir quantification
    cd quantification
    STUDY=$(basename $(dirname "$PWD"))
    echo "(10/10) Gene count"
    find ../aligned -name "*sortedByReadname.out.bam" | while read file ; do xabse=$(basename $file); echo "echo -e \"$file\t\$(samtools view -h -b -q 10 -f 2 "$file" | samtools view -c)\" >> ${STUDY}.totalreadcount.txt" >> ${STUDY}.totalreadcount.commands.txt ; done ;
    parallel_GNU -j $(wc -l ${STUDY}.totalreadcount.commands.txt | awk '{print $1}') < ${STUDY}.totalreadcount.commands.txt
    featureCounts_input_files=`ls ../aligned/*sortedByReadname.out.bam`
    featureCounts -a /bar/yliang/genomes/private/gencode.v36.primary_assembly.annotation.gtf -o ${STUDY}_featurecounts.txt -t exon -T 8 -s 1 -g gene_name -p -B -C $featureCounts_input_files 2> ${STUDY}_featurecounts.log
    python3 /bar/yliang/tricks/rna_pipe_featureCounts_to_z-score.py --input ${STUDY}_featurecounts.txt --size ${STUDY}.totalreadcount.txt
fi






