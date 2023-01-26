#!/bin/bash -e

# programmer: Holden
# This is a pipline to do the initial steps with wgbs data up to assembly

# Arguments to bash script to decide on various parameters
MAXJOBS=6 #Cores per CPU
MAXJOBSALIGN=4 # Cores for each alignment
FASTQFOLDER="./fastq"
READ1EXTENSION="R1.fastq.gz"

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

cd $FASTQFOLDER/..

echo "parallel job limit: ${MAXJOBS}"
echo "parallel job limit for alignment: ${MAXJOBSALIGN}"
echo "fastqfile location: ${FASTQFOLDER}"
echo "read 1 extension: ${READ1EXTENSION}"
echo "reference genome for bismark: ${GENOMEDIR}"

#ml cutadapt STAR FastQC stringtie python2 multiqc picard
ml bismark bwa samtools FastQC python2 multiqc preseq picard

## save softwares_version information
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

## start processing
mkdir trimmed
cd trimmed
#find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o Trimed_"$xbase" -p Trimed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 
find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "gunzip -c "$file" | fastx_trimmer -f 10 -l 150 -Q33 -z -o Trimmed_"$xbase" " >> 1_cutadaptcommands.txt ; done;
find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "gunzip -c "${file/R1/R2}" | fastx_trimmer -f 16 -l 150 -Q33 -z -o Trimmed_"${xbase/R1/R2}" " >> 1_cutadaptcommands.txt ; done;
echo "(1/5) Trimming Reads"
#parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt


find ./ -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc "$file"" >> 1b_fastqc_commands.txt; done ;
echo "(1b/5) FastQC run"
#parallel_GNU -j $MAXJOBS < 1b_fastqc_commands.txt

ml python2
multiqc .
ml python3

cd ..
mkdir lambda_aligned
cd lambda_aligned
#bismark --bowtie2 --multicore 8 -N 1 --score_min L,0,-0.2 -o ./lambda_aligned /scratch/twlab/hlee/genomes/lambda -1 ./trimmedreads/fastx_trim/${part1}'_R1_fastxtrimmed.fq.gz' -2 ./trimmedreads/fastx_trim/${part1}'_R2_fastxtrimmed.fq.gz'
find ../trimmed -maxdepth 1 -name "Trimmed_*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "bismark --bowtie2 --multicore "$MAXJOBSALIGN" -N 1 --score_min L,0,-0.2 -o ./ --genome "$LAMBDAGENOMEDIR" -1 "$file" -2 "${file/R1/R2}" &> ./"$xbase"_bismark_align.log" >> 2_lambda_alignCommands.txt ; done;
echo "(2a/5) Aligning Reads to lambda genome"
#parallel_GNU -j $MAXJOBS < 2_lambda_alignCommands.txt


cd ..
mkdir aligned
cd aligned
#bismark -I $MIN_INSERT -X $MAX_INSERT --parallel $CPU -o $outdir --temp_dir ./tmp --gzip --nucleotide_coverage $genome_dir -1 $inputR1 -2 $inputR2 &>/scratch/twlab/jjang/TE_antigen/WGBS/aligned/${part1}'_bismark.log'
find ../trimmed -name "Trimmed_*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "bismark --bowtie2 -I $MIN_INSERT -X $MAX_INSERT --multicore "$MAXJOBSALIGN" -o ./ --temp_dir ./"$xbase"_tmp --gzip --nucleotide_coverage --genome "$GENOMEDIR" -1 "$file" -2 "${file/R1/R2}" &> ./"$xbase"_bismark_align.log " >> 2_alignCommands.txt ; done ;
echo "(2b/5) Aligning Reads to human genome"
#parallel_GNU -j $MAXJOBS < 2_alignCommands.txt

#samtools sort Trimed_SCC90_BRep1_WGBS_R1_bismark_bt2_pe.bam -o Trimed_SCC90_BRep1_WGBS_R1_bismark_bt2_pe.sorted.bam -O bam -@ 2 -T Trimed_SCC90_BRep1_WGBS_R1_bismark_bt2_pe
find ./ -name "*bam" | while read file; do xbase=$(basename $file ".bam"); echo "samtools sort "$file" -o "$xbase".sorted.bam -O bam -@ 5 -T "$xbase""  >> 3_sortingCommands.txt ; done ;
echo "(3/5) Sorting bam files"
#parallel_GNU -j $MAXJOBS < 3_sortingCommands.txt

#deduplicate_bismark -p --output_dir $outdir --bam $bam &>$log_dedup
find ./ -name "*sorted.bam" | while read file; do xbase=$(basename $file); echo "deduplicate_bismark -p --bam "$file" &> ./"$xbase"_bismark_deduplication.log "  >> 4_deduplicationCommands.txt ; done ;
echo "(4/5) Deduplication"
#parallel_GNU -j $MAXJOBS < 4_deduplicationCommands.txt

cd ..
mkdir methylcall
cd methylcall
#bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report -o $outdir2 --gzip --parallel $CPU $bam_dedup &>$log_methx
find ../aligned -name "*.deduplicated.bam" | while read file; do xbase=$(basename $file); echo "bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report -o ../methylcall --gzip --multicore "$MAXJOBSALIGN" "$file" &> ../methylcall/"$xbase"_bismark_methylcall.log " >> 5_extractorCommands.txt ; done ;
echo "(5/5) Methylation Extraction"
#parallel_GNU -j $MAXJOBS < 5_extractorCommands.txt

#cd ..
#mkdir complexity_QC
#cd complexity_QC
#find ../aligned -name "*sorted.bam" | while read file; do xbase=$(basename $file ".sorted.bam"); echo "preseq c_curve -pe -bam -o "$xbase".c_curve.txt "$file" &> ./"$xbase".preseq_c_curve.log" >> 6a_preseq_c_curveCommands.txt ; done ;
#echo "(6a/) preseq_c_curve"
#parallel_GNU -j $MAXJOBS < 6a_preseq_c_curveCommands.txt
#python3 /bar/yliang/tricks/preseq_c_curve_plot.py --input $PWD --output $(basename $(dirname "$PWD")).c_curve.png
#
#find ../aligned -name "*sorted.bam" | while read file; do xbase=$(basename $file ".sorted.bam"); echo "preseq lc_extrap -pe -bam -o "$xbase".lc_extrap.txt "$file" &> ./"$xbase".preseq_lc_extrap.log" >> 6b_preseq_lc_extrapCommands.txt ; done ;
#echo "(6b/) preseq_lc_extrap"
#parallel_GNU -j $MAXJOBS < 6b_preseq_lc_extrapCommands.txt
#python3 /bar/yliang/tricks/preseq_lc_extrap_plot.py --input $PWD --output $(basename $(dirname "$PWD")).lc_extrap.png
#
#find ../aligned -name "*sorted.bam" | while read file; do xbase=$(basename $file ".sorted.bam"); echo "preseq lc_extrap -pe -bam -o "$xbase".lc_extrapdefects.txt "$file" -defects &> ./"$xbase".preseq_lc_extrapdefects.log" >> 6c_preseq_lc_extrapdefectsCommands.txt ; done ;
#echo "(6b/) preseq_lc_extrapdefects"
#parallel_GNU -j $MAXJOBS < 6c_preseq_lc_extrapdefectsCommands.txt
#python3 /bar/yliang/tricks/preseq_lc_extrapdefects_plot.py --input $PWD --output $(basename $(dirname "$PWD")).lc_extrapdefects.png
#
#find ../aligned -name "*sorted.bam" | while read file ; do xbase=$(basename $file ".sorted.bam");echo "java -jar \$PICARD EstimateLibraryComplexity I="$file" O="$xbase".picard_estimateComplexity.txt &> ./"$xbase".picard_estimateComplexity.log" >> 6d_picard_estimateComplexityCommands.txt ; done ;
#echo "(6c/) picard_estimateComplexity"
#parallel_GNU -j $MAXJOBS < 6d_picard_estimateComplexityCommands.txt
#grep -A 1 ESTIMATED_LIBRARY_SIZE *picard_estimateComplexity.txt > $(basename $(dirname "$PWD")).picard_estimateComplexity.txt

#ls *bam | while read file ; do samtools index $file ; samtools idxstats $file > ${file}_idxstats ; done ;
#
#cd ..
#mkdir rseqc
#cd rseqc
#
#find ../aligned -name "*bam" | while read file ; do xbase=$(basename $file) ; echo "geneBody_coverage.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38.HouseKeepingGenes.bed" >> 3_rseqQCcommands.txt ; echo "read_distribution.py -i "$file" -r ~/reference/rseqc/hg38_Gencode_V28.bed > "${xbase%.*}".readdistribution.txt" >> 3_rseqQCcommands.txt ; echo "junction_saturation.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38_Gencode_V28.bed" >> 3_rseqQCcommands.txt ; echo "bam2wig.py -s "~/reference/mapability/hg38.chrom.sizes" -i "$file" -o "${xbase%.*}" -u" >> 3_rseqQCcommands.txt; done ;
#
#echo "(3/5) Quality Control and Visualization"
#parallel_GNU -j $MAXJOBS < 3_rseqQCcommands.txt
#
#find ../aligned -name "*bam" | while read file ; do echo "java -jar $PICARD EstimateLibraryComplexity I="$file" O="${file%.*}"_duplication_stats.txt" >> estimateComplexity.txt ; done ;
#
#echo "(3.5/5) Estimate Complexity"
#parallel_GNU -j $MAXJOBS < estimateComplexity.txt
#
#cd ..
#mkdir assembled
#cd assembled
#
#find ../aligned -name "*bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -q 255 -h "$file" | stringtie - -o "${xbase%.*}".gtf -p 4 -m 100 -c 1 --fr" >> 4_assembleCommands.txt ; done ;
#
#echo "(4/5) Quality Control and Visualization"
#parallel_GNU -j $MAXJOBS < 4_assembleCommands.txt
#
#cd ..
#echo "(5/5) multiQC"
#
#multiqc .
#
#