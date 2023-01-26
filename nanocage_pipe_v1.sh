#!/bin/bash
#!/usr/bin/awk
#####################################################################
# This is a pipeline for nanoCAGE-seq alignment 
# programmer: Holden Liang (adapted from Juheon's pipeline)
# /scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/pipe_align_nanoCAGE.sh
#####################################################################
# Arguments to bash script to decide on various parameters
MAXJOBS=6 # Cores per CPU
MAXJOBS_STAR=2 # STAR requires a lot of resources so it is not revommended to have more than 2 alignment running at the same time.
CAPFILTER=0.15 # default threshold for G-Cap filtering
FASTQFOLDER="/scratch/yliang/HNSCC/data/nanocage_keratinocyte_QC/fastq"
READ1EXTENSION="R1.fastq.gz"

ADAPTOR="xiaoyun" ## adaptor sequence, if have "ACACAG" in all reads, choose xiaoyun, otherwise, choose juheon
GENOME="hg38" ## hg38 or mm10 

#####################################################################
REF="/bar/genomes/hg38/hg38.fa"
#REF_STAR="/bar/genomes/hg38/STAR"
REF_STAR="/bar/yliang/genomes/private/STAR_index_hg38_gencodeV22"
GTF_FILE="/bar/yliang/genomes/private/gencode.v29.annotation.gtf"
FGBIO="/bar/yliang/softwares/fgbio/fgbio/target/scala-2.13/fgbio-1.4.0-29132fb-SNAPSHOT.jar"
TagdustOutput2UMI="/bar/yliang/tricks/nanocage_pipe_TagdustOutput2UMI.py"
QC1="/bar/yliang/tricks/nanocage_pipe_QC1.R"

PRIMER_5END_RC=CCCTATA            #Reverse-complement of TSO: TATAGGG
PRIMER_3END_RC=ACATCTCCGAGCCCACGAGAC  #Reverse-complement of transposon: GTCTCGTGGGCTCGGAGATGT

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
    -a|--adaptor)
    ADAPTOR="$2"
    shift # past argument
    ;;
    -g|--genome)
    GENOME="$2"
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

ml cutadapt STAR FastQC 

########################################################################
## 1-1. Trmming of UMI and TSO at 5' end of R1
cd $FASTQFOLDER/..
mkdir trimmed
cd trimmed

# xiaoyun ADAPTOR: tagdust -t $MAXJOBS_STAR -1 O:N -2 S:ACACAG -3 F:NNNNNNNNN -4 S:TATAGGG -5 R:N -show_finger_seq "$file" "${file/R1/R2}" -o .
# juheon ADAPTOR: tagdust -t $MAXJOBS_STAR -1 O:N -2 F:NNNNNNNNN -3 S:TATAGGG -4 R:N -show_finger_seq "$file" "${file/R1/R2}" -o .
if [[ $ADAPTOR == "xiaoyun" ]]; then
    find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "tagdust -t $MAXJOBS_STAR -1 O:N -2 S:ACACAG -3 F:NNNNNNNNN -4 S:TATAGGG -5 R:N -show_finger_seq "$file" "${file/R1/R2}" -o TD_"$xbase" &> "$xbase"_tagdust.log" >> 1a_tagdust_commands.txt ; done ; 
elif [[ $ADAPTOR == "juheon" ]]; then
    find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "tagdust -t $MAXJOBS_STAR -1 O:N -2 F:NNNNNNNNN -3 S:TATAGGG -4 R:N -show_finger_seq "$file" "${file/R1/R2}" -o TD_"$xbase" &> "$xbase"_tagdust.log" >> 1a_tagdust_commands.txt ; done ; 
fi
echo "(1a/) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1a_tagdust_commands.txt

########################################################################
## 1-2. get UMI
# python2 /bar/yliang/tricks/nanocage_pipe_TagdustOutput2UMI.py 
find ../trimmed -maxdepth 1 -name "TD_*${READ1EXTENSION}_READ1.fq" | while read file ; do xbase=$(basename $file) ; echo "python2 $TagdustOutput2UMI $file "$file"_UMI" >> 1b_UMI_commands.txt; done;
echo "(1b/) Getting UMI"
parallel_GNU -j $MAXJOBS < 1b_UMI_commands.txt

gzip *fq

########################################################################
## 1-3. Trimming of adapter sequences
## Used arguments: -a: adapter at the 3' end of READ 1, -A: adapter at the 3' end of READ 2, -m: minimum length of reads after the adapter trimming
## If the fragment is short, Read 1 may contain 3' end seqeuncing primer and vice versa. This step is to get rid of such primer sequences
## Should I use smaller number for the minimum length? For the READ 1, read size after UMI and TS oligo removal is 59 bp. If it contains PCR_REV_RC, it will be likely discarded
## https://cutadapt.readthedocs.io/en/stable/recipes.html, "Trimming (amplicon-) primers from both ends of paired-end reads"
## cutadapt does not do reverse complement on its own

# cutadapt -j $MAXJOBS_STAR -a "$PRIMER_3END_RC;min_overlap=10" -A "$PRIMER_5END_RC;min_overlap=7" -m 50 -o Trimed_"$xbase" -p Trimed_"${xbase/R1/R2}" "$file" "${file/R1/R2}"
find ../trimmed -maxdepth 1 -name "TD_*${READ1EXTENSION}_READ1.fq.gz" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a \"$PRIMER_3END_RC;min_overlap=10\" -A \"$PRIMER_5END_RC;min_overlap=7\" --minimum-length 50 -o Trimmed_"$xbase" -p Trimmed_"${xbase/READ1/READ2}" "$file" "${file/READ1/READ2}" > "$xbase"_cutadapt.log" >> 1c_cutadapt_commands.txt ; done ; 
echo "(1c/) Trimming Reads"
parallel_GNU -j $MAXJOBS < 1c_cutadapt_commands.txt

########################################################################
## 1-4. FastQC
find ../trimmed -maxdepth 1 -name "Trimmed*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc $file" >> 1d_fastqc_commands.txt; done ;
echo "(1d/) FastQC run"
parallel_GNU -j $MAXJOBS < 1d_fastqc_commands.txt

########################################################################
## 2. Alignment using STAR
# --twopassMode Basic: Better to find unannottated splice junctions
# --chimOutType WithinBAM SoftClip: Chimeric alignments will be included in the output BAM 
# --sjdbGTFfile $GTFFILE: Known splice junctions are helpful for the alignment to splice junctions
cd ..
mkdir aligned
cd aligned

find ../trimmed -name "Trimmed*_READ1.fq.gz" | while read file; do xbase=$(basename $file); echo "STAR --runThreadN 4 --runMode alignReads --twopassMode Basic --chimOutType WithinBAM SoftClip --sjdbGTFfile $GTF_FILE --genomeDir $REF_STAR --outFileNamePrefix "${xbase%.*}" --readFilesIn "$file" "${file/READ1/READ2}" --readFilesCommand zcat --outSAMtype BAM  SortedByCoordinate" >> 2_align_commands.txt ; done ;
echo "(2/) Aligning Reads with STAR"
parallel_GNU -j $MAXJOBS_STAR < 2_align_commands.txt

########################################################################
## 3. Mark UMIs
# java -jar $FGBIO AnnotateBamWithUmis -i ./Trimmed_${xbase/_UMI/Aligned.sortedByCoord.out.bam} -f $file -o ./Trimmed_${xbase/_UMI/Aligned.sortedByCoord.out.UMI.bam}
find ../trimmed -name "*UMI" | while read file ; do xbase=$(basename $file) ; echo "java -jar $FGBIO AnnotateBamWithUmis -i ./Trimmed_${xbase/_UMI/Aligned.sortedByCoord.out.bam} -f $file -o ./Trimmed_${xbase/_UMI/Aligned.sortedByCoord.out.UMI.bam}" >> 3_markUMI_commands.txt; done ;
echo "(3/) Mark UMIs"
parallel_GNU -j $MAXJOBS < 3_markUMI_commands.txt

########################################################################
## 4. Deduplication using UMIs
ml python3
# umi_tools dedup --paired --extract-umi-method=tag --umi-tag=RX -I $file -L $xbase".log" -E $xbase".err" --output-stats=$xbase".stats" -S ${xbase/bam/dedup.bam}
find ../aligned -name "*UMI.bam" | while read file ; do xbase=$(basename $file) ; echo "umi_tools dedup --paired --extract-umi-method=tag --umi-tag=RX -I $file -L ${xbase/bam/dedup.bam}".log" -E ${xbase/bam/dedup.bam}".err" --output-stats=${xbase/bam/dedup.bam}".stats" -S ${xbase/bam/dedup.bam}" >> 4_deduplication_commands.txt; done ;
echo "(4/) Deduplication"
parallel_GNU -j $MAXJOBS < 4_deduplication_commands.txt

########################################################################
## 5. Deduplication using UMIs
# samtools view -h -f 0x2 -F 0x900 -q 255 -o $DEDUP_FILTBAM $DEDUPBAM
find ../aligned -name "*dedup.bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -h -f 0x2 -F 0x900 -q 255 -b -o ${xbase/bam/filtered.bam} $file" >> 5_filter_commands.txt; done ;
echo "(5/) Filter reads"
parallel_GNU -j $MAXJOBS < 5_filter_commands.txt

########################################################################
## 6. rseqc
cd ..
mkdir rseqc
cd rseqc

find ../aligned -name "*filtered.bam" | while read file ; do xbase=$(basename $file) ; echo "geneBody_coverage.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38.HouseKeepingGenes.bed" >> 6_rseqQCcommands.txt ; echo "read_distribution.py -i "$file" -r ~/reference/rseqc/hg38_Gencode_V28.bed > "${xbase%.*}".readdistribution.txt" >> 6_rseqQCcommands.txt ; echo "junction_saturation.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38_Gencode_V28.bed" >> 6_rseqQCcommands.txt ; echo "bam2wig.py -s "~/reference/mapability/hg38.chrom.sizes" -i "$file" -o "${xbase%.*}" -u" >> 6_rseqQCcommands.txt; done ;


########################################################################
## 7. multiQC
cd ..
echo "(6/) multiQC"
/bar/yliang/anaconda3/bin/multiqc .

########################################################################
## 8. QC1--check reads and peaks distribution
mkdir QC
cd QC

ml R/3.6.1
# Rscript $QC1 -wd $PWD -input $file -samplename $xbase -genome $GENOME 
find ../aligned -name "*filtered.bam" | while read file ; do xbase=$(basename $file "_READ1.fqAligned.sortedByCoord.out.UMI.dedup.filtered.bam") ; echo "Rscript $QC1 -wd $PWD -input $file -samplename $xbase -genome $GENOME &> $xbase.QC1.log" >> 8_QC1_commands.txt; done ;
echo "(8/) Filter reads"
parallel_GNU -j $MAXJOBS < 8_QC1_commands.txt

########################################################################
## 9. peak calling
## adapted from /scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/pipe_bam_to_peak.sh
cd ../aligned

# bamtoCTSS: /scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/pipe_bam_to_peak.sh
# python2 /bar/yliang/tricks/nanocage_bam2CTSS.py $filtered.bam .
find . -name "*filtered.bam" | while read file; do xbase=$(basename $file); echo "python2 /bar/yliang/tricks/nanocage_bam2CTSS.py $file ." >> 9a_PeakCalling_commands.txt; done;
echo "(9/) Peak Calling"
parallel_GNU -j $MAXJOBS < 9a_PeakCalling_commands.txt

find . -name "*bam.CTSS" | while read file; do echo "sort -k1,1 -V -k2,2n $file | egrep -v chrUn | egrep -v chrM |  egrep -v _random | egrep -v _alt > ${file/CTSS/nochrM.CTSS}" >> 9b_PeakCalling_commands.txt; done;
find . -name "*_unannotatedG.CTSS" | while read file; do echo "sort -k1,1 -V -k2,2n $file | egrep -v chrUn | egrep -v chrM |  egrep -v _random | egrep -v _alt > ${file/CTSS/nochrM.CTSS}" >> 9b_PeakCalling_commands.txt; done;
parallel_GNU -j $MAXJOBS < 9b_PeakCalling_commands.txt

find . -name "*nochrM.CTSS" | while read file ; do echo "awk -v OFS=\"\t\" '{print \$1,\$2,\$2+1,\".\",\$4,\$3}' "$file" | sort -k1,1V -k2,2n > "$file".bdg ; bgzip "$file".bdg ; tabix -p bed "$file".bdg.gz" >> 9c_PeakCalling_commands.txt ; done ;
parallel_GNU -j $MAXJOBS < 9c_PeakCalling_commands.txt

STUDY=$(basename $(dirname "$PWD"))
ml R/3.6.1
# peak calling: /scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/CAGEr_nakulpeakcalling_nrPassThreshold_paraclu.R
echo "Rscript /bar/yliang/tricks/nanocage_CAGEr_nakulpeakcalling_nrPassThreshold_paraclu.R . $STUDY.CTSS" >> 9d_PeakCalling_commands.txt
parallel_GNU -j 1 < 9d_PeakCalling_commands.txt

# post-processing (Gcap filter): /scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/post_process_CAGEpeak_v2.R
echo "Rscript /bar/yliang/tricks/nanocage_post_process_CAGEpeak_v2.R $STUDY.CTSS $CAPFILTER $STUDY " >> 9e_PeakCalling_commands.txt
parallel_GNU -j 1 < 9e_PeakCalling_commands.txt


################################################################################################################################################################################################################################################################################################
##!/bin/bash
## programmer: Nakul
## This is a pipline to do the initial steps with CAGE-seq data up to assembly and then the main pipeline that we have can be run
#
## Arguments to bash script to decide on various parameters
#MAXJOBS=6 #Cores per CPU
#MAXJOBSSTAR=2 #STAR requires a lot of resources so it is not revommended to have more than 2 alignment running at the same time. 
#FASTQFOLDER="./"
#READ1EXTENSION="R1.fastq.gz"
#
#while [[ $# > 1 ]]
#do
#key="$1"
#case $key in
#    -j|--maxjobs)
#    MAXJOBS="$2"
#    shift # past argument
#    ;;
#    -s|--maxstarjobs)
#    MAXJOBSSTAR="$2"
#    shift # past argument
#    ;;
#    -f|--fastqfolder)
#    FASTQFOLDER="$2"
#    shift # past argument
#    ;;
#    -r|--read1extension)
#    READ1EXTENSION="$2"
#    shift # past argument
#    ;;
#    --default)
#    DEFAULT=YES
#    ;;
#    *)
#            # unknown option
#    ;;
#esac
#shift # past argument or value
#done
#
#echo "parallel job limit: ${MAXJOBS}"
#echo "parallel job limit STAR: ${MAXJOBSSTAR}"
#echo "fastqfile location: ${FASTQFOLDER}"
#echo "read 1 extension: ${READ1EXTENSION}"
#
#ml FastQC fastx-toolkit samtools STAR cutadapt multiqc
#
#mkdir removeforwardtag
#cd removeforwardtag
#
#find $FASTQFOLDER -maxdepth 1 -name  "*${READ1EXTENSION}"  | while read file ; do xbase=$(basename $file); echo "gunzip -c "$file" | fastx_trimmer -Q33 -f 17 -z > "$xbase >> 1_trim18commands.txt ; done ;
#
#echo "(1/8) Remove Forward Tag Reads and Link Reverse"
#parallel_GNU -j $MAXJOBS < 1_trim18commands.txt
#
#find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION/R1/R2}" | while read file ; do ln -s $file ; done ;
#
#cd ..
#
#mkdir trimmed
#cd trimmed
#
#find ../removeforwardtag -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A TGCTGGAGACCTTCAGTTCGACTAGTGTAG --minimum-length 50  -o Trimed_"$xbase" -p Trimed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 2_cutadaptcommands.txt ; done ; 
#
#echo "(2/8) Trimming Reads"
#parallel_GNU -j $MAXJOBS < 2_cutadaptcommands.txt
#
#find ../trimmed -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc "$file >> 3_fastqc_commands.txt; done ;
#
#echo "(3/8) FastQC run"
#parallel_GNU -j $MAXJOBS < 3_fastqc_commands.txt
#
#cd ..
#mkdir aligned
#cd aligned
#
#find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "STAR  --runMode alignReads  --runThreadN 4  --genomeDir ~/reference/STAR_index_hg38_gencodeV22/ --readFilesIn "$file" "${file/R1/R2}" --readFilesCommand zcat --outFileNamePrefix "${xbase%.*}" --outSAMtype BAM   SortedByCoordinate   --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMheaderHD @HD VN:1.4 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33  --alignIntronMax 500000  --alignMatesGapMax 1000000 --twopassMode Basic" >> 4_alignCommands.txt ; done ;
#
#echo "(4/8) Aligning Reads with STAR"
#parallel_GNU -j $MAXJOBSSTAR < 4_alignCommands.txt
#
#ls *bam | while read file ; do samtools index $file ; samtools idxstats $file > ${file}_idxstats ; done ;
#
#cd ..
#mkdir rseqc
#cd rseqc
#
#find ../aligned -name "*bam" | while read file ; do xbase=$(basename $file) ; echo "geneBody_coverage.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38.HouseKeepingGenes.bed" >> 5_rseqQCcommands.txt ; echo "read_distribution.py -i "$file" -r ~/reference/rseqc/hg38_Gencode_V28.bed > "${xbase%.*}".readdistribution.txt" >> 5_rseqQCcommands.txt ; echo "junction_saturation.py -i "$file" -o "${xbase%.*}" -r ~/reference/rseqc/hg38_Gencode_V28.bed" >> 5_rseqQCcommands.txt ; echo "bam2wig.py -s "~/reference/mapability/hg38.chrom.sizes" -i "$file" -o "${xbase%.*}" -u" >> 5_rseqQCcommands.txt; done ;
#
#echo "(5/8) Quality Control and Visualization"
#parallel_GNU -j $MAXJOBS < 5_rseqQCcommands.txt
#
#cd ..
#mkdir nsorted
#cd nsorted
#
#find ../aligned -name "*bam" | while read file ; do xbase=$(basename $file) ; echo "samtools sort -n -o "${xbase/sortedByCoord/sortedByName}" "$file >> 6_sortedNameCommands.txt ; done ;
#
#echo "(6/8) Quality Control and Visualization"
#parallel_GNU -j $MAXJOBS < 6_sortedNameCommands.txt
#
#find . -name "*bam" | while read file ; do echo "bam2CTSS.py "$file >> 7_createCTSSCommands.txt ; done ;
#
#echo "(7/8) CTSS creation and visualization"
#parallel_GNU -j $MAXJOBS < 7_createCTSSCommands.txt
#
#find . -name "*.CTSS" | while read file ; do echo "CTSStoBedgraph.py "$file" | sort -k1,1V -k2,2n > "$file".bdg ; bgzip "$file".bdg ; tabix -p bed "$file".bdg.gz" >> 7a_ctssVisualizeCommands.txt ; done ;
#
#cd ..
#echo "(8/8) multiQC"
#
#multiqc .

