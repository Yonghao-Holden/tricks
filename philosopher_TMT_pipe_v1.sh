#!/bin/bash
# programmer: Holden
#####################################################################
# This is a pipline to perform initial analysis on processed MS spectrum data (mzML files) using philosopher pipeline
# https://github.com/Nesvilab/philosopher/wiki/Pipeline-mode-for-TMT-analysis
#####################################################################
## step 1. download data from CPTAC and organized in a specific format
## https://github.com/esacinc/PDC-Public/tree/master/tools/downloadPDCData
## python3 downloadPDCData.py <PDC File manifest.csv>
python3 /bar/yliang/tricks/downloadPDCData.py PDC_file_manifest_11032022_114341.csv

## step 2. move files from our server to htcf for further processing
rsync -avFP 

philosopher workspace --clean




## step 3. set up 
philosopher pipeline --config params/philosopher.yml 01CPTAC_HNSCC_W_JHU_20190507 02CPTAC_HNSCC_W_JHU_20190513 2>> PDC000221_HNSCC_mini.err &>> PDC000221_HNSCC_mini.log






find ../aligned -name "*sortedByCoord.out.bam" | while read file ; do xbase=$(basename $file) ; echo "samtools view -q 255 -h "$file" | stringtie - -o "${xbase%.*}".gtf -p 4 -m 100 -c 1 --fr" >> 8_assembleCommands.txt ; done ;


# https://github.com/Nesvilab/philosopher/issues/135
philosopher labelquant --bestpsm --brand itraq --dir . --level 2 --plex 4 --removelow