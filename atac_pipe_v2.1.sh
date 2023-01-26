#!/bin/bash
#!/usr/bin/awk
# programmer: Nakul & Holden
# Run it under the folder `<ATAC folder>/`
#####################################################################
# This is a full ATAC-seq pipeline to get all the way to peaks from a set of files.
# IDR and other calculations can be done at another time 
#####################################################################
# default settings
MAXJOBS=5 #Cores per CPU
MAXJOBS_PICRAD=2 #Cores per CPU
FASTQFOLDER="/scratch/yliang/Epitherapy_3D/data/atac-seq_AllSample/fastq"
READ1EXTENSION="R1.fastq.gz"
MACSQVALUE=".05"
myQC="yes"

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

cd $FASTQFOLDER/..

ml cutadapt bowtie2 FastQC samtools picard macs2 multiqc python2 R/3.4.1

## save softwares_version information
echo "cutadapt version:" >> softwares_version.txt
cutadapt --version >> softwares_version.txt
echo "bowtie2 version:" >> softwares_version.txt
bowtie2 --version >> softwares_version.txt
echo "FastQC version:" >> softwares_version.txt
fastqc --version >> softwares_version.txt
echo "samtools version:" >> softwares_version.txt
samtools --version >> softwares_version.txt
echo "picard_CollectInsertSizeMetrics version:" >> softwares_version.txt
java -jar $PICARD CollectInsertSizeMetrics --version 2>> softwares_version.txt
echo "macs2 version:" >> softwares_version.txt
macs2 --version 2>> softwares_version.txt
echo "multiQC version:" >> softwares_version.txt
multiqc --version >> softwares_version.txt
echo "python2 version:" >> softwares_version.txt
python2 --version 2>> softwares_version.txt
echo "R version: 3.4.1" >> softwares_version.txt

mkdir trimmed
cd trimmed


find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --quality-cutoff=15,10 --minimum-length=36 -o Trimmed_"$xbase" -p Trimmed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 
echo "(1/17) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt 2> 1_cutadaptcommands.err &> 1_cutadaptcommands.log

find ../trimmed -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc $file" >> fastqc_commands.txt; done ;
echo "(1b/17) FastQC run"
parallel_GNU -j $MAXJOBS < fastqc_commands.txt

cd ..
mkdir aligned
cd aligned

find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "bowtie2 --very-sensitive -p 4 -X 2000 -x /bar/yliang/genomes/private/hg38_bowtie2_index/hg38_bw -1 "$file" -2 "${file/R1/R2}" | samtools view -u - | samtools sort - > "${xbase%.*}".bam ; samtools index "${xbase%.*}".bam" >> 2_alignCommands.txt ; done ;
echo "(2/17) Aligning Reads with bowtie2"
parallel_GNU -j $MAXJOBS < 2_alignCommands.txt 2> 2_alignCommands.err &> 2_alignCommands.log

find . -name "*.bam" | while read file ; do samtools idxstats $file > ${file}_idxstats ; done ;

## Remove chrM and other chromosomes
find . -name "*fastq.bam" | while read file ; do echo "samtools view -b "$file" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort > "${file%.*}"_nochrM.bam ; samtools index "${file%.*}"_nochrM.bam" >> 3_remove_chrM.txt ; done ;
echo "(3/17) chrM Removal"
parallel_GNU -j $MAXJOBS < 3_remove_chrM.txt 2> 3_remove_chrM.err &> 3_remove_chrM.log

## Remove duplicates and get stats
## this step requires a lot of memory --> only 3 jobs at a time 
find . -name "*_nochrM.bam" | while read file ; do echo "java -jar $PICARD MarkDuplicates I="$file" O="${file%.*}"_nodup.bam M="$file"_dups.txt REMOVE_DUPLICATES=true" >> 4_removeDuplicates.txt ; done ;
echo "(4/17) Duplicate Marking and Removal"
parallel_GNU -j 3 < 4_removeDuplicates.txt 2> 4_removeDuplicates.err &> 4_removeDuplicates.log

##Calculate flagstat
find . -name "*_nodup.bam" | while read file ; do echo "samtools flagstat "$file" > "$file"_flag.stats" >> 5_flagstat_calc.txt ; done ;
echo "(5/17) Calculating Flag Statistics"
parallel_GNU -j $MAXJOBS < 5_flagstat_calc.txt 2> 5_flagstat_calc.err &> 5_flagstat_calc.log

##Extract Properly Paired and Uniquely Mapped
find . -name "*_nodup.bam" | while read file ; do echo "samtools view -h -b -q 10 -f 2 "$file" > "${file%.*}"_qffilter.bam" >> 6_extractProperlyPaired.txt ; done ;
echo "(6/17) Removing Unpaired and Low Quality Reads"
parallel_GNU -j $MAXJOBS < 6_extractProperlyPaired.txt 2> 6_extractProperlyPaired.err &> 6_extractProperlyPaired.log

rm *fastq.bam
rm *nochrM.bam
rm *nodup.bam
rm *bai

##Index the most recent file
find . -name "*_qffilter.bam" | while read file ; do echo "samtools index "$file >> 7_indexCommands.txt ; done ;
echo "(7/17) Indexing Files"
parallel_GNU -j $MAXJOBS < 7_indexCommands.txt 2> 7_indexCommands.err &> 7_indexCommands.log

##Calculate insert size distribution
find . -name "*_qffilter.bam" | while read file ; do echo "java -jar $PICARD CollectInsertSizeMetrics I="$file" O="${file%.*}"_insertsize.txt H="${file%.*}"_insertsize.pdf" >> 8_insertSizeDistribution.txt ; done ;
echo "(8/17) Calculating Insert Size Distribution"
parallel_GNU -j $MAXJOBS < 8_insertSizeDistribution.txt 2> 8_insertSizeDistribution.err &> 8_insertSizeDistribution.log

##Generate Tag aligned files with shifted coordinates
find . -name "*qffilter.bam" | while read file ; do echo "samtools sort -n "$file" | bedtools bamtobed -bedpe -i stdin | awk -v OFS=\"\t\" '{if(\$9==\"+\"){print \$1,\$2+4,\$6+4}else if(\$9==\"-\"){print \$1,\$2-5,\$6-5}}' > "${file%.*}".tagAlign" >> 9_tagAlignGenerate.txt ; done ;
echo "(9/17) Shifting BAM coordinates to tagAlign"
parallel_GNU -j $MAXJOBS < 9_tagAlignGenerate.txt 2> 9_tagAlignGenerate.err &> 9_tagAlignGenerate.log

##Call peaks
find . -name "*tagAlign" | while read file ; do echo "macs2 callpeak -t "$file" -f BED -n "${file%.*}" -q "$MACSQVALUE" --nomodel --shift 37 --extsize 73 -B --keep-dup all --call-summits" >> 10_macs2Call.txt ; done ;
echo "(10/17) MACS2 Peak Calling"
parallel_GNU -j $MAXJOBS < 10_macs2Call.txt 2> 10_macs2Call.err &> 10_macs2Call.log

find . -name "*tagAlign" | while read file ; do echo "macs2 callpeak -t "$file" -f BED -n "${file%.*}"_uniqpeak -q "$MACSQVALUE" --nomodel --shift 37 --extsize 73 -B --keep-dup all" >> 10_macs2Call_uniqpeak.txt ; done ;
echo "(10/17) MACS2 Peak Calling"
parallel_GNU -j $MAXJOBS < 10_macs2Call_uniqpeak.txt 2> 10_macs2Call_uniqpeak.err &> 10_macs2Call_uniqpeak.log

##Remove blacklisted sequences
#This is from ENCODE and is only about 38 sequences
find . -name "*.narrowPeak" | while read file ; do echo "bedtools intersect -v -a "$file" -b /bar/yliang/genomes/private/hg38.blacklist.bed > "${file%.*}"_noBL.narrowPeak" >> 11_removeBlacklist.txt ; done ;
echo "(11/17) Remove Blacklist Regions"
parallel_GNU -j $MAXJOBS < 11_removeBlacklist.txt 2> 11_removeBlacklist.err &> 11_removeBlacklist.log

##Generate bigWig to be used for browser visualization
find . -name "*_treat_pileup.bdg" | while read file ; do echo "sort -k1,1 -k2,2n "$file" > "${file%.*}"_sorted.bdg ; bedGraphToBigWig "${file%.*}"_sorted.bdg /bar/yliang/genomes/private/hg38_all_chromosomes.size "${file%.*}"_sorted.bw" >> 12_generateBW.txt ; done ;
echo "(12/17) Generate BigWig Files"
parallel_GNU -j $MAXJOBS < 12_generateBW.txt 2> 12_generateBW.err &> 12_generateBW.log


##Generate sorted and indexed bed file that can be used on WashU browser
find . -name "*_peaks.narrowPeak" | while read file ; do echo "awk '{print \$1\"\t\"\$2\"\t\"\$3}' $file | sort -k1,1 -k2,2n > "${file}"_sorted ; bedtools merge -i "${file}"_sorted > "${file}"_merge_sorted ; bgzip "${file}"_merge_sorted ; tabix -p bed "${file}"_merge_sorted.gz" >> 13_processBed.txt ; done ;
echo "(13/17) Generate peak bed files"
parallel_GNU -j $MAXJOBS <  13_processBed.txt 2> 13_processBed.err &> 13_processBed.log

##########################################################################################
## Some more QC analysis
find . -name "*_peaks.narrowPeak" | while read file ; do echo "/bar/yliang/myapps/bin/ataqv --ignore-read-groups --peak-file "$file" --tss-file /bar/yliang/genomes/private/hg38.tss.refseq.bed human "${file/_peaks.narrowPeak/.bam} >> 14_atacQCcommands.txt ; done ;
echo "(14/17) ATACQV Commands"
parallel_GNU -j $MAXJOBS <  14_atacQCcommands.txt 2> 14_atacQCcommands.err &> 14_atacQCcommands.log

echo "(15/17) ATACQV Visualization Commands"
FILES=$(ls *.ataqv.json | tr '\n' ' ')
mkarv webapp $FILES

echo "(16/17) Fragment Length Distribution"
find . -name "*qffilter.bam" | while read file ; do echo "samtools view "$file" | awk '\$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > "$file"_fragdist; Rscript /bar/yliang/softwares/Rscripts/fragSizeDist.R "$file"_fragdist" >> 16_fragDistCommands.txt ; done ;
parallel_GNU -j $MAXJOBS <  16_fragDistCommands.txt 2> 16_fragDistCommands.err &> 16_fragDistCommands.log

echo "(17/17) MultiQC Commands"
cd ..
/bar/yliang/anaconda3/bin/multiqc .

mv multiqc_report.html $(basename "$PWD").multiqc_report.html

echo "Pipeline started at"$START" and ended at "$END
echo ""
echo "Now you can view the QC using the MultiQC reports and also the ataqv visualizations that have been generated (copy the webapp folder to a public link)."
echo ""
echo "In addition, the bigwig and bed files can be visualized on the WashU Epigenome Browser. Using the following commands to make a json hub:"
echo "makedatahub.py .bw https://wangftp.wustl.edu/~nshah/ucsf/TE_Antigen_Project/Pilot3/ATAC/ Trimed_WangT_ATAC_ Pilot3 bigWig green ATAC > ATAChub_Pilot2_bw.json"
echo "makedatahub.py _sorted.gz https://wangftp.wustl.edu/~nshah/ucsf/TE_Antigen_Project/Pilot3/ATAC/ Trimed_WangT_ATAC_ Pilot3 bed black Peaks > ATAChub_Pilot2_bed.json"


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

