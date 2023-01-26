#!/bin/bash
#!/usr/bin/awk
# programmer: Holden Liang
# Run it under the folder `centrifuge/`
#####################################################################
# This is a pipeline to centrifuge sequence reads to different genomes.
# https://github.com/DaehwanKimLab/centrifuge
#####################################################################
# default settings
MAXJOBS=6 #Cores per CPU
TRIMMEDFOLDER="../trimmed/"
READ1EXTENSION="R1.fastq.gz"
GENOME="/bar/yliang/softwares/centrifuge/centrifuge/indices/p+h+v"
END="PE"

# Arguments to bash script to decide on various parameters
while [[ $# > 1 ]]
do
key="$1"
case $key in
    -j|--maxjobs)
    MAXJOBS="$2"
    shift # past argument
    ;;
    -t|--trimmedfolder)
    TRIMMEDFOLDER="$2"
    shift # past argument
    ;;
    -r|--read1extension)
    READ1EXTENSION="$2"
    shift # past argument
    ;;
    -g|--refgenome)
    GENOME="$2"
    shift # past argument
    ;;
    -e|--end)
    END="$2"
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
echo "fastq file location: ${TRIMMEDFOLDER}"
echo "read1extension: ${READ1EXTENSION}"
echo "ref_genome: ${GENOME}"

#module load matlab python3

cd $TRIMMEDFOLDER/..

mkdir centrifuge
cd centrifuge

## save softwares_version information
echo "Centrifuge version:" >> softwares_version.txt
centrifuge --version >> softwares_version.txt

## Centrifuge trimmed reads
#centrifuge -q -p 8 -x /bar/nshah/reference/centrifuge/p+h+v -1 /scratch/nakul/HNSCCollab/20200310_celllines/trimmed/Trimed_UMSCC-47_HPV_R1.fastq.gz -2 /scratch/nakul/HNSCCollab/20200310_celllines/trimmed/Trimed_UMSCC-47_HPV_R2.fastq.gz -S Trimed_UMSCC-47_HPV_R1.fastq.gz_centrifuge_classification --report-file Trimed_UMSCC-47_HPV_R1.fastq.gz_centrifuge_report.txt
if [ ${END} = "PE" ]
then
    find ${TRIMMEDFOLDER} -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file ${READ1EXTENSION}) ; echo "centrifuge -q -p 8 -x ${GENOME} -1 $file -2 ${file/R1/R2} -S "$xbase"centrifuge_classification --report-file "$xbase"centrifuge_report.txt"  >> 1_centrifuge_commands.txt ; done ;
elif [ ${END} = "SE" ]
then
    READ1EXTENSION=".fastq.gz"
    find ${TRIMMEDFOLDER} -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=$(basename $file ${READ1EXTENSION}) ; echo "centrifuge -q -p 8 -x ${GENOME} -U $file -S "$xbase"_centrifuge_classification --report-file "$xbase"_centrifuge_report.txt"  >> 1_centrifuge_commands.txt ; done ;
else
    echo "input PE or SE for --end flag"
fi

echo "centrifuge"
parallel_GNU -j $MAXJOBS < 1_centrifuge_commands.txt


