#!/bin/bash
#!/usr/bin/awk
# programmer: Holden
#####################################################################
# This is a pipline to perform initial analysis on RNA-seq data including de novo assemble and gene expression quantification, which will be ready to use for TEProf2 and DESeq2 directly, respectively.
# run it like this: bash /bar/yliang/tricks/rna_pipe.sh -f <path_to_fastq, like /scratch/yliang/rna/fastq>
#####################################################################
# Update log:
# 1. pangenome library prep kit is rf. add an option for pangenome data (rf data)
#####################################################################
# default settings
MAXJOBS=6 #Cores per CPU
MAXJOBSSTAR=2 #STAR requires a lot of resources so it is not revommended to have more than 2 alignment running at the same time. 
FASTQFOLDER="./" ## please use absolute path to run 
READ1EXTENSION="R1.fastq.gz"
QUANTIFICATION="yes"
TWOPASS="no" ## run two pass alignment for STAR
STRAND="rf" ## old kit (universal) we use is fr. But in most cases like pangenome data, encode data and data generated using NEB kit, it's rf

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
    -t|--twopass)
    TWOPASS="$2"
    shift # past argument
    ;;
    -st|--strand)
    STRAND="$2"
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

if [ -d "./trimmed" ] 
then
    rm -rf trimmed aligned rseqc assembled quantification softwares_version.txt
fi

ml cutadapt STAR FastQC python2 picard samtools

## save softwares_version information
echo "cutadapt version:" >> softwares_version.txt
cutadapt --version >> softwares_version.txt
echo "STAR version:" >> softwares_version.txt
STAR --version >> softwares_version.txt
echo "FastQC version:" >> softwares_version.txt
fastqc --version >> softwares_version.txt
echo "stringtie version:" >> softwares_version.txt
/bar/yliang/myapps/bin/stringtie --version >> softwares_version.txt
echo "multiQC version:" >> softwares_version.txt
/bar/yliang/anaconda3/bin/multiqc --version >> softwares_version.txt
echo "picard_EstimateLibraryComplexity version:" >> softwares_version.txt
java -jar $PICARD EstimateLibraryComplexity --version 2>> softwares_version.txt
echo "samtools version:" >> softwares_version.txt
samtools --version >> softwares_version.txt
echo "featureCounts version:" >> softwares_version.txt
/bar/yliang/myapps/bin/featureCounts -v 2>> softwares_version.txt


cd $FASTQFOLDER

find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o ${xbase%.*}_fastqc "$file >> 2_fastqc_commands.txt; done ;
echo "(2a/10) FastQC run on the raw fastq file"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 2_fastqc_commands.txt

cd ..
mkdir trimmed
cd trimmed

find $FASTQFOLDER -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file) ; echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o Trimmed_"$xbase" -p Trimmed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt ; done ; 
echo "(1/10) Trimming Reads"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 1_cutadaptcommands.txt

find ../trimmed -maxdepth 1 -name "*.gz"  | while read file ; do xbase=$(basename $file) ; mkdir ${xbase%.*}_fastqc ; echo "fastqc -o "${xbase%.*}"_fastqc "$file >> 2_fastqc_commands.txt; done ;
echo "(2b/10) FastQC run on the trimmed fastq file"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 2_fastqc_commands.txt

cd ..
mkdir aligned
cd aligned

if [[ $TWOPASS == "no" ]]
then
    find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "STAR  --runMode alignReads  --runThreadN 4  --genomeDir /bar/yliang/genomes/private/STAR_index_hg38_gencodeV36/ --readFilesIn "$file" "${file/R1/R2}" --readFilesCommand zcat --outFileNamePrefix "${xbase%.*}" --outSAMtype BAM   SortedByCoordinate   --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMheaderHD @HD VN:1.4 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33  --alignIntronMax 500000  --alignMatesGapMax 1000000 --twopassMode Basic" >> 3_alignCommands.txt ; done ;
    echo "(3/10) Aligning Reads with STAR"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSSTAR < 3_alignCommands.txt
fi

## this session follow the gdc guideline as i handle CPTAC RNA-seq data deposited on GDC
## https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
## they used gencode v36
if [[ $TWOPASS != "no" ]]
then
    find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); mkdir ${xbase/R1.fastq.gz/1stpass_output}; echo "STAR --genomeDir /bar/yliang/genomes/private/STAR_index_hg38_gencodeV36/ --readFilesIn "$file" "${file/R1/R2}" --runThreadN 4 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None --outFileNamePrefix ./${xbase/R1.fastq.gz/1stpass_output}/ 2> ${xbase/R1.fastq.gz/1stpass_output}/1st_pass.err &> ${xbase/R1.fastq.gz/1stpass_output}/1st_pass.log" >> 3_alignCommands_1stpass.txt ; done ;
    echo "(3/10) Aligning Reas with STAR (first pass)"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSSTAR < 3_alignCommands_1stpass.txt

    find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "STAR --runMode genomeGenerate --genomeDir ${xbase/R1.fastq.gz/1stpass_output} --genomeFastaFiles /bar/yliang/genomes/private/GRCh38.primary_assembly.genome.fa --sjdbOverhang 100 --runThreadN 4 --sjdbFileChrStartEnd ${xbase/R1.fastq.gz/1stpass_output}/SJ.out.tab --outTmpDir ${xbase/R1.fastq.gz/genome_index_temp} 2> ${xbase/R1.fastq.gz/1stpass_output}/genome_indexing.err &> ${xbase/R1.fastq.gz/1stpass_output}/genome_indexing.log "  >> 3_alignCommands_generate_index.txt ; done ;
    echo "(3/10) Aligning Reas with STAR (generating index)"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSSTAR < 3_alignCommands_generate_index.txt

    find ../trimmed -name "*${READ1EXTENSION}" | while read file; do xbase=$(basename $file); echo "STAR --runMode alignReads --genomeDir ${xbase/R1.fastq.gz/1stpass_output} --readFilesIn "$file" "${file/R1/R2}" --runThreadN 4 --outFileNamePrefix "${xbase%.*}" --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 2> ${xbase/R1.fastq.gz/2nd_pass.err} &> ${xbase/R1.fastq.gz/2nd_pass.log}" >> 3_alignCommands_2ndpass.txt;done;
    echo "(3/10) Aligning Reas with STAR (second pass)"
    /bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSSTAR < 3_alignCommands_2ndpass.txt
fi

find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "samtools index $file ; samtools idxstats $file > ${file}_idxstats" >> 4_indexingCommands.txt; done;
echo "(4/10) index bam files"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 4_indexingCommands.txt

#find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "samtools index $file ; samtools idxstats $file > ${file}_idxstats ; samtools sort $file -O bam -@ 3 -n -T ${file}_temp > ${file/sortedByCoord/sortedByReadname} " >> 4_sortingCommands.txt; done;
#echo "(4/10) Sort bam files by readname"
#/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 4_sortingCommands.txt

if [[ $STRAND == "fr" ]]; then
    find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "/bar/yliang/anaconda3/bin/bamCoverage -b $file -o ${file/bam/reverse.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand forward --scaleFactor -1 2> 5_bam2bigwigCommand.log" >> 5_bam2bigwigCommand.txt; done;
    find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "/bar/yliang/anaconda3/bin/bamCoverage -b $file -o ${file/bam/forward.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand reverse 2> 5_bam2bigwigCommand.log" >> 5_bam2bigwigCommand.txt; done;
elif [[ $STRAND == "rf" ]] ; then
    find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "/bar/yliang/anaconda3/bin/bamCoverage -b $file -o ${file/bam/forward.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand forward 2> 5_bam2bigwigCommand.log" >> 5_bam2bigwigCommand.txt; done;
    find ../aligned -name "*sortedByCoord.out.bam" | while read file; do xbase=$(basename $file); echo "/bar/yliang/anaconda3/bin/bamCoverage -b $file -o ${file/bam/reverse.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand reverse --scaleFactor -1 2> 5_bam2bigwigCommand.log" >> 5_bam2bigwigCommand.txt; done;
fi
echo "(5/10) Convert bam to bigwig"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 5_bam2bigwigCommand.txt

cd ..
mkdir rseqc
cd rseqc

find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "geneBody_coverage.py -i "$file" -o "${xbase%.*}" -r /bar/yliang/genomes/private/rseqc/hg38.HouseKeepingGenes.bed" >> 6_rseqQCcommands.txt ; echo "read_distribution.py -i "$file" -r /bar/yliang/genomes/private/rseqc/hg38_Gencode_V28.bed > "${xbase%.*}".readdistribution.txt" >> 6_rseqQCcommands.txt ; echo "junction_saturation.py -i "$file" -o "${xbase%.*}" -r /bar/yliang/genomes/private/rseqc/hg38_Gencode_V28.bed" >> 6_rseqQCcommands.txt >> 6_rseqQCcommands.txt; done ;
echo "(6/10) Quality Control and Visualization"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 6_rseqQCcommands.txt

find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "java -jar \$PICARD EstimateLibraryComplexity I="$file" O="${file%.*}"_duplication_stats.txt" >> 7_estimateComplexity.txt ; done ;
echo "(7/10) Estimate Complexity"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBSSTAR < 7_estimateComplexity.txt

cd ..
mkdir assembled
cd assembled

if [[ $STRAND == "fr" ]]; then
    find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -q 255 -h "$file" | /bar/yliang/myapps/bin/stringtie - -o "${xbase%.*}".gtf -p 4 -m 100 -c 1 --fr" >> 8_assembleCommands.txt ; done ;
elif [[ $STRAND == "rf" ]] ; then
    find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -q 255 -h "$file" | /bar/yliang/myapps/bin/stringtie - -o "${xbase%.*}".gtf -p 4 -m 100 -c 1 --rf" >> 8_assembleCommands.txt ; done ;
fi
echo "(8/10) Quality Control and Visualization"
/bar/yliang/myapps/bin/parallel_GNU -j $MAXJOBS < 8_assembleCommands.txt

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
    find ../aligned -name "*.sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -q 255 -h "$file" | /bar/yliang/myapps/bin/stringtie - -o "${xbase%.*}".gtf -e -b "${xbase%.*}".stats -p 4 -m 100 -c 1 -G /bar/yliang/genomes/private/gencode.v36.primary_assembly.annotation.gtf" >> quantification_commands.txt ; done ;
    parallel_GNU -j $(wc -l quantification_commands.txt | awk '{print $1}') < quantification_commands.txt
    python3 /bar/yliang/tricks/compile_stringtie_quantification.py --output ${STUDY}_stringtie_quantification.txt
fi


