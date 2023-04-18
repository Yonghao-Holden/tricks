#!/bin/bash -e

# programmer: Holden
# This is a pipline to do the initial steps with wgbs data up to assembly
# run it like this: bash /bar/yliang/tricks/wgbs_pipe.sh -f <path_to_fastq, like /scratch/yliang/wgbs/fastq>
########################################################################
# Arguments to bash script to decide on various parameters
MAXJOBS=4 #Cores per CPU
MAXJOBSALIGN=4 # Cores for each alignment
FASTQFOLDER="./fastq" ## please use absolute path to run 
READ1EXTENSION="R1.fastq.gz"
COMPLEXITY="YES"
########################################################################
GENOMEDIR="/bar/yliang/genomes/public/hg38/bismark"
LAMBDAGENOMEDIR="/bar/yliang/genomes/private/lambda_bismark"
MIN_INSERT=0
MAX_INSERT=2000
QC1="/bar/yliang/tricks/wgbs_pipe_QC1.R"

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
    -c|--complexity)
    COMPLEXITY="$2"
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

ml bismark bwa samtools FastQC picard

cd $FASTQFOLDER/..
STUDY=$(basename $PWD)
## save softwares_version information
echo "fastx_toolkit version:" >> softwares_version.txt
/bar/yliang/myapps/bin/fastx_trimmer -h >> softwares_version.txt
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
find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "gunzip -c "$file" | /bar/yliang/myapps/bin/fastx_trimmer -f 10 -l 150 -Q33 -z -o Trimmed_"$xbase" " >> 1_cutadaptcommands.txt ; done;
find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "gunzip -c "${file/R1/R2}" | /bar/yliang/myapps/bin/fastx_trimmer -f 16 -l 150 -Q33 -z -o Trimmed_"${xbase/R1/R2}" " >> 1_cutadaptcommands.txt ; done;
echo "(1/5) Trimming Reads"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt

## 1b. fastqc
find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; mkdir Trimmed_${xbase%.*}_fastqc ; echo "fastqc -o Trimmed_${xbase%.*}_fastqc "$file >> 1b_fastqc_commands.txt; done ;
echo "(1b/5) FastQC run"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 1b_fastqc_commands.txt

########################################################################
## 2a. align to lambda genome
cd ..
mkdir lambda_aligned
cd lambda_aligned
#bismark --bowtie2 --multicore 8 -N 1 --score_min L,0,-0.2 -o ./lambda_aligned /scratch/twlab/hlee/genomes/lambda -1 ./trimmedreads/fastx_trim/${part1}'_R1_fastxtrimmed.fq.gz' -2 ./trimmedreads/fastx_trim/${part1}'_R2_fastxtrimmed.fq.gz'
find ../trimmed -maxdepth 1 -name "Trimmed_*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "bismark --bowtie2 --multicore "$MAXJOBSALIGN" -N 1 --score_min L,0,-0.2 -o ./ --genome "$LAMBDAGENOMEDIR" -1 "$file" -2 "${file/R1/R2}" &> ./"$xbase"_bismark_align.log" >> 2a_lambda_alignCommands.txt ; done;
echo "(2a/5) Aligning Reads to lambda genome"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 2a_lambda_alignCommands.txt
grep -E "Mapping efficiency:" *PE_report.txt >> lambda_alignment_report.txt
grep -E "C methylated in" *PE_report.txt >> lambda_alignment_report.txt

########################################################################
## 2b. align to human genome
cd ..
mkdir aligned
cd aligned
#bismark -I $MIN_INSERT -X $MAX_INSERT --parallel $CPU -o $outdir --temp_dir ./tmp --gzip --nucleotide_coverage $genome_dir -1 $inputR1 -2 $inputR2 &>/scratch/twlab/jjang/TE_antigen/WGBS/aligned/${part1}'_bismark.log'
find ../trimmed -name "Trimmed_*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "bismark --bowtie2 -I $MIN_INSERT -X $MAX_INSERT --multicore "$MAXJOBSALIGN" -o ./ --temp_dir ./"$xbase"_tmp --gzip --nucleotide_coverage --genome "$GENOMEDIR" -1 "$file" -2 "${file/R1/R2}" &> ./"$xbase"_bismark_align.log " >> 2b_alignCommands.txt ; done ;
echo "(2b/5) Aligning Reads to human genome"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 2b_alignCommands.txt

########################################################################
## 3. deduplication
#deduplicate_bismark -p --output_dir $outdir --bam $bam &>$log_dedup
find ./ -name "*pe.bam" | while read file; do xbase=$(basename $file); echo "deduplicate_bismark -p --bam "$file" &> ./"$xbase"_bismark_deduplication.log "  >> 3_deduplicationCommands.txt ; done ;
echo "(3/5) Deduplication"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 3_deduplicationCommands.txt

########################################################################
## 4. call methylation
#bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report -o $outdir2 --gzip --parallel $CPU $bam_dedup &>$log_methx
find ./ -name "*.deduplicated.bam" | while read file; do xbase=$(basename $file); echo "bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report -o ./ --gzip --bedGraph --multicore "$MAXJOBSALIGN" "$file" &> ./"$xbase"_bismark_methylcall.log " >> 4_extractorCommands.txt ; done ;
echo "(4/5) Methylation Extraction"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 4_extractorCommands.txt

########################################################################
## 5-8. complexity check
if [[ $COMPLEXITY != "no" ]]
then
    find ./ -name "*pe.bam" | while read file; do xbase=$(basename $file ".bam"); echo "samtools sort "$file" -o "$xbase".sorted.bam -O bam -@ 5 -T "$xbase"" >> 6_sortingCommands.txt; done ;
    echo "(6/) sorting bam files"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 6_sortingCommands.txt
    cd ..
    mkdir complexity_QC
    cd complexity_QC
    find ../aligned -name "*sorted.bam" | while read file; do xbase=$(basename $file ".sorted.bam"); echo "/bar/yliang/myapps/bin/preseq c_curve -pe -bam -o "$xbase".c_curve.txt "$file" &> ./"$xbase".preseq_c_curve.log" >> 7_preseq_c_curveCommands.txt ; done ;
    echo "(7/) preseq_c_curve"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 7_preseq_c_curveCommands.txt
    python3 /bar/yliang/tricks/preseq_c_curve_plot.py --input $PWD --output ${STUDY}.c_curve.png
    find ../aligned -name "*sorted.bam" | while read file; do xbase=$(basename $file ".sorted.bam"); echo "/bar/yliang/myapps/bin/preseq lc_extrap -pe -bam -o "$xbase".lc_extrap.txt "$file" &> ./"$xbase".preseq_lc_extrap.log" >> 8a_preseq_lc_extrapCommands.txt ; done ;
    echo "(8a/) preseq_lc_extrap"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 8a_preseq_lc_extrapCommands.txt
    python3 /bar/yliang/tricks/preseq_lc_extrap_plot.py --input $PWD --output ${STUDY}.lc_extrap.png --defect no
    find ../aligned -name "*sorted.bam" | while read file; do xbase=$(basename $file ".sorted.bam"); echo "/bar/yliang/myapps/bin/preseq lc_extrap -pe -bam -o "$xbase".lc_extrapdefects.txt "$file" -defects &> ./"$xbase".preseq_lc_extrapdefects.log" >> 8b_preseq_lc_extrapdefectsCommands.txt ; done ;
    echo "(8b/) preseq_lc_extrapdefects"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 8b_preseq_lc_extrapdefectsCommands.txt
    python3 /bar/yliang/tricks/preseq_lc_extrap_plot.py --input $PWD --output ${STUDY}.lc_extrapdefects.png --defect yes
    find ../aligned -name "*sorted.bam" | while read file ; do xbase=$(basename $file ".sorted.bam");echo "java -jar \$PICARD EstimateLibraryComplexity I="$file" O="$xbase".picard_estimateComplexity.txt &> ./"$xbase".picard_estimateComplexity.log" >> 9_picard_estimateComplexityCommands.txt ; done ;
    echo "(9/) picard_estimateComplexity"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSALIGN < 9_picard_estimateComplexityCommands.txt # don't run too many samples simultaneously, picard takes up a lot of memory
    for file in *picard_estimateComplexity.txt
    do
        xbase=$(basename $file ".R1_bismark_bt2_pe.picard_estimateComplexity.txt")
        number=`grep -A 1 ESTIMATED_LIBRARY_SIZE $file | tail -n 1`
        echo -e "$xbase\t$number" >> ${STUDY}.picard_estimateComplexity.txt
    done
    python3 /bar/yliang/tricks/picard_complexity_plot.py --input $STUDY.picard_estimateComplexity.txt --output $STUDY.picard_complexity.png
fi

########################################################################
## 9. report gereration
bismark2report
bismark2summary
cd ..
/bar/yliang/anaconda3/bin/multiqc .
mv multiqc_report.html $(basename "$PWD").multiqc_report.html

########################################################################
## 10. sample correlation
#cd ../aligned
#SAMPLENUM=`ll *.deduplicated.bismark.cov.gz|wc -l`
#zcat *.deduplicated.bismark.cov.gz | cut -f 1-3 | sort -Vk1,2 | uniq -c | awk -v OFS="\t" -v SAMPLENUM=$SAMPLENUM '{if($1==SAMPLENUM) print $2,$3,$4}' | grep -P "chr1\t|chr2\t|chr3\t|chr4\t|chr5\t|chr6\t|chr7\t|chr8\t|chr9\t|chr10\t|chr11\t|chr12\t|chr13\t|chr14\t|chr15\t|chr16\t|chr17\t|chr18\t|chr19\t|chr20\t|chr21\t|chr22\t" | sort -Vk1,2 |uniq > ${STUDY}_common_autosome_CpG_across_samples.bed
#
#find ../aligned -name "*.deduplicated.bismark.cov.gz" | while read file; do echo "zcat $file | bedtools intersect -a stdin -b ${STUDY}_common_autosome_CpG_across_samples.bed -u -f 1| sort -Vk1,2 | uniq > ${file/.gz/.commonCpG.bed}" >> 10_get_commonCpG_info_commands.txt; done;
#echo "(10/) get common CpG information"
#parallel_GNU -j $MAXJOBS < 10_get_commonCpG_info_commands.txt
#
#Rscript $QC1 --wd $PWD --input _R1_bismark_bt2_pe.deduplicated.bismark.cov.commonCpG.bed --output $STUDY
#
#gzip *.commonCpG.bed
#rm ${STUDY}_common_autosome_CpG_across_samples.bed



























