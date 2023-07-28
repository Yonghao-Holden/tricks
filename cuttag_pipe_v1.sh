#!/bin/bash
#!/usr/bin/awk
# programmer: Holden
#####################################################################
# This is a full CUT&TAG pipeline to get all the way to peaks from a set of files.
# Setting of Bowtie2 and MACS2 comes from here: https://www.nature.com/articles/s41467-019-09982-5#Sec8 
#####################################################################
# bash /bar/yliang/tricks/chip_pipe.sh -f /scratch/yliang/share/for_xiaoyun/GSE124557_CUTTAG_K562_H3K27Ac_H3K4Me3/fastq -mq yes -g hg38 -i CUTTAG_K562_input_R1.fastq.gz -ml 15
# default settings
MAXJOBS=6 #Cores per CPU
MAXJOBS_PICRAD=2 #Cores per CPU
FASTQFOLDER="/scratch/yliang/HNSCC/data/atac-seq_HIOEC/fastq/"
READ1EXTENSION="R1.fastq.gz"
MACSQVALUE=".05"
myQC="yes"
GENOME="hg38"
INPUT="input_N3_R1.fastq.gz" # input fastq file name
MIN_LENGTH=36

# Arguments to bash script to decide on various parameters
while [[ $# > 1 ]]
do
key="$1"
case $key in
    -j|--maxjobs)
    MAXJOBS="$2"
    shift # past argument
    ;;
    -jp|--maxjobspicard)
    MAXJOBS_PICRAD="$2"
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
    -m|--macs2qthresh)
    MACSQVALUE="$2"
    shift # past argument
    ;;
    -mq|--myqc)
    myQC="$2"
    shift # past argument
    ;;
    -g|--genome)
    GENOME="$2"
    shift # past argument
    ;;
    -i|--input)
    INPUT="$2"
    shift # past argument
    ;;
    -ml|--minlength)
    MIN_LENGTH="$2"
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
echo "fastqfile location: ${FASTQFOLDER}"
echo "read 1 extension: ${READ1EXTENSION}"
echo "macs2 q-value: ${MACSQVALUE}"

echo "Make sure that bowtie2, fastqc, samtools, picard, macs2, multiqc, ataaqv, fragSizeDist.R, and R (ggplot2, grid, gridExtra) have all been loaded in the environment."

module load cutadapt bowtie2 FastQC samtools picard macs2 multiqc python2 R/3.4.1

cd $FASTQFOLDER/..
mkdir trimmed
cd trimmed


find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --quality-cutoff=15,10 --minimum-length=$MIN_LENGTH -o Trimmed_"$xbase" -p Trimmed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 

echo "(1/17) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt

find ../trimmed -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc $file" >> fastqc_commands.txt; done ;

echo "(1b/17) FastQC run"
parallel_GNU -j $MAXJOBS < fastqc_commands.txt

cd ..
mkdir aligned
cd aligned

if [[ $GENOME == "hg38" ]]; then
    find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "bowtie2 --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x /bar/yliang/genomes/private/hg38_bowtie2_index/hg38_bw -1 "$file" -2 "${file/R1/R2}" | samtools view -u - | samtools sort - > "${xbase%.*}".bam ; samtools index "${xbase%.*}".bam" >> 2_alignCommands.txt ; done ;
elif [[ $GENOME == "mm10" ]]; then
    find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "bowtie2 --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x /bar/yliang/genomes/private/mm10_bowtie2_index/mm10 -1 "$file" -2 "${file/R1/R2}" | samtools view -u - | samtools sort - > "${xbase%.*}".bam ; samtools index "${xbase%.*}".bam" >> 2_alignCommands.txt ; done ;
fi

echo "(2/17) Aligning Reads with bowtie2"
parallel_GNU -j $MAXJOBS < 2_alignCommands.txt

find . -name "*.bam" | while read file ; do samtools idxstats $file > ${file}_idxstats ; done ;

##Remove chrM and other chromosomes
find . -name "*fastq.bam" | while read file ; do echo "samtools view -b "$file" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort > "${file%.*}"_nochrM.bam ; samtools index "${file%.*}"_nochrM.bam" >> 3_remove_chrM.txt ; done ;

echo "(3/17) chrM Removal"
parallel_GNU -j $MAXJOBS < 3_remove_chrM.txt

##Remove duplicates and get stats

find . -name "*_nochrM.bam" | while read file ; do echo "java -jar $PICARD MarkDuplicates I="$file" O="${file%.*}"_nodup.bam M="$file"_dups.txt REMOVE_DUPLICATES=true" >> 4_removeDuplicates.txt ; done ;

echo "(4/17) Duplicate Marking and Removal"
parallel_GNU -j $MAXJOBS < 4_removeDuplicates.txt

##Calculate flagstat
find . -name "*_nodup.bam" | while read file ; do echo "samtools flagstat "$file" > "$file"_flag.stats" >> 5_flagstat_calc.txt ; done ;

echo "(5/17) Calculating Flag Statistics"
parallel_GNU -j $MAXJOBS < 5_flagstat_calc.txt

##Extract Properly Paired and Uniquely Mapped
find . -name "*_nodup.bam" | while read file ; do echo "samtools view -h -b -q 10 -f 2 "$file" > "${file%.*}"_qffilter.bam" >> 6_extractProperlyPaired.txt ; done ;

echo "(6/17) Removing Unpaired and Low Quality Reads"
parallel_GNU -j $MAXJOBS < 6_extractProperlyPaired.txt

rm *fastq.bam
rm *nochrM.bam
rm *nodup.bam
rm *bai

##Index the most recent file
find . -name "*_qffilter.bam" | while read file ; do echo "samtools index "$file >> 7_indexCommands.txt ; done ;

echo "(7/17) Indexing Files"
parallel_GNU -j $MAXJOBS < 7_indexCommands.txt

##Calculate insert size distribution
find . -name "*_qffilter.bam" | while read file ; do echo "java -jar $PICARD CollectInsertSizeMetrics I="$file" O="${file%.*}"_insertsize.txt H="${file%.*}"_insertsize.pdf" >> 8_insertSizeDistribution.txt ; done ;

echo "(8/17) Calculating Insert Size Distribution"
parallel_GNU -j $MAXJOBS < 8_insertSizeDistribution.txt

##Call peaks
find . -name "*_qffilter.bam" | while read file ; do if [[ $file != "./Trimmed_${INPUT/.gz/_nochrM_nodup_qffilter.bam}" ]]; then echo "macs2 callpeak -t "$file" -c Trimmed_"${INPUT/.gz/_nochrM_nodup_qffilter.bam}" -f BAMPE -n "${file%.*}" -p 0.00001 --keep-dup all --call-summits" >> 10_macs2Call.txt ;fi ; done ;

echo "(10/17) MACS2 Peak Calling"
parallel_GNU -j $MAXJOBS < 10_macs2Call.txt

find . -name "*_qffilter.bam" | while read file ; do if [[ $file != "./Trimmed_${INPUT/.gz/_nochrM_nodup_qffilter.bam}" ]]; then echo "macs2 callpeak -t "$file" -c Trimmed_"${INPUT/.gz/_nochrM_nodup_qffilter.bam}" -f BAMPE -n "${file%.*}"_uniqpeak -p 0.00001 --keep-dup all" >> 10_macs2Call_uniqpeak.txt; fi ; done ;

echo "(10/17) MACS2 Peak Calling"
parallel_GNU -j $MAXJOBS < 10_macs2Call_uniqpeak.txt 

##Remove blacklisted sequences
#This is from ENCODE and is only about 38 sequences (hg38)
##########################################################################################
if [[ $GENOME == "hg38" ]]; then
    find . -name "*.narrowPeak" | while read file ; do echo "bedtools intersect -v -a "$file" -b /bar/yliang/genomes/private/hg38.blacklist.bed > "${file%.*}"_noBL.narrowPeak" >> 11_removeBlacklist.txt ; done ;
elif [[ $GENOME == "mm10" ]]; then
    find . -name "*.narrowPeak" | while read file ; do echo "bedtools intersect -v -a "$file" -b /bar/yliang/genomes/private/mm10.blacklist.bed > "${file%.*}"_noBL.narrowPeak" >> 11_removeBlacklist.txt ; done ;
fi

echo "(11/17) Remove Blacklist Regions"
parallel_GNU -j $MAXJOBS < 11_removeBlacklist.txt

##Generate bigWig to be used for browser visualization
find . -name "*_nodup_qffilter.bam" | while read file ; do echo "/bar/yliang/anaconda3/bin/bamCoverage --bam $file -o ${file/bam/bw} -of bigwig --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads" >> 12_generateBW.txt ; done ; # http://genomewiki.ucsc.edu/index.php/Mm10_Genome_size_statistics non-N bases

echo "(12/17) Generate BigWig Files"
parallel_GNU -j $MAXJOBS < 12_generateBW.txt
##########################################################################################
##Generate sorted and indexed bed file that can be used on WashU browser
find . -name "*noBL.narrowPeak" | while read file ; do echo "awk '{print \$1\"\t\"\$2\"\t\"\$3}' $file | sort -k1,1 -k2,2n > "${file}"_sorted ; bedtools merge -i "${file}"_sorted > "${file}"_merge_sorted ; bgzip "${file}"_merge_sorted ; tabix -p bed "${file}"_merge_sorted.gz" >> 13_processBed.txt ; done ;

echo "(13/17) Generate peak bed files"
parallel_GNU -j $MAXJOBS <  13_processBed.txt 

##########################################################################################
## Some more QC analysis
#find . -name "*_peaks.narrowPeak" | while read file ; do echo "ataqv --ignore-read-groups --peak-file "$file" --tss-file /bar/yliang/genomes/private/hg38.tss.refseq.bed human "${file/_peaks.narrowPeak/.bam} >> 14_atacQCcommands.txt ; done ;
#
#echo "(14/17) ATACQV Commands"
#parallel_GNU -j $MAXJOBS <  14_atacQCcommands.txt

#echo "(15/17) ATACQV Visualization Commands"
#FILES=$(ls *.ataqv.json | tr '\n' ' ')
#mkarv webapp $FILES

echo "(16/17) Fragment Length Distribution"
find . -name "*qffilter.bam" | while read file ; do echo "samtools view "$file" | awk '\$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > "$file"_fragdist; Rscript /bar/yliang/softwares/Rscripts/fragSizeDist.R "$file"_fragdist" >> 16_fragDistCommands.txt ; done ;
parallel_GNU -j $MAXJOBS <  16_fragDistCommands.txt

echo "(17/17) MultiQC Commands"
cd ..
/bar/yliang/anaconda3/bin/multiqc .

if [[ $myQC == "yes" ]]
then
    mkdir myQC
    cd myQC
    STUDY=$(basename $(dirname "$PWD"))
    echo "(18/18)myQC"
    ## 1. percentage of read under peaks
    echo "percentage of read under peaks"
    find ../aligned -name "*_qffilter.bam" | while read file; do peak_file=${file/.bam/_uniqpeak_peaks_noBL.narrowPeak}; output_file=${file/.bam/.bamintersectpeak} ; output_file_base=$(basename $output_file) ;echo "bedtools intersect -a $peak_file -b $file -bed -loj -wa -wb > $output_file_base" >> 18c_RUP_1.commands.txt; done
    parallel_GNU -j 1 < 18c_RUP_1.commands.txt
    find . -name "*bamintersectpeak" | while read file; do xbase=$(basename $file); echo "echo -e \"${xbase/.bamintersectpeak/.bam}\t\$(cut -f 14 $file | sort | uniq -c | wc -l)\" >> ${STUDY}.RUPcount.txt " >> 18c_RUP_2.commands.txt; done
    parallel_GNU -j $MAXJOBS < 18c_RUP_2.commands.txt
    find ../aligned -name "*_qffilter.bam" | while read file; do xbase=$(basename $file); echo "echo -e \"${xbase}\t\$(samtools view -c $file)\" >> ${STUDY}.totalreadcount.txt" >> 18c_RUP_3.commands.txt; done
    #find ../aligned -name "*_qffilter.bam" | while read file; do xbase=$(basename $file); echo "samtools view -c $file > ${xbase/.bam/.readcount.txt}" >> 18c_RUP_2.commandds.txt; done
    parallel_GNU -j $MAXJOBS < 18c_RUP_3.commands.txt
    rm *bamintersectpeak

    ## 2. read location distribution
    echo "(19/19) read location distribution"

    ## 3. peak location distribution
    echo "(20/20) peak location distribution"
    find ../aligned -name "*_uniqpeak_peaks_noBL.narrowPeak" | while read file; do xbase=$(basename $file) ;echo "bedtools intersect -a $file -b /bar/yliang/genomes/private/GENCODE_V36/GENCODE_v36_HL_genic_annotation.bed -f 0.8 -loj -wa -wb > ${xbase/_uniqpeak_peaks_noBL.narrowPeak/.PeakGenicIntersection}" >> 20a_commands.txt; done

    find . -name "*PeakGenicIntersection" | while read file; do echo $file; cut -f 14 $file | sort | uniq -c ; done

    ## 4. peak size distribution
fi



#sendmail holdenleungyh@gmail.com < /bar/yliang/tricks/mail_ATAC.txt
