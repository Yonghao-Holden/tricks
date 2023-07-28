#!/bin/bash
# programmer: Holden Liang
# this script uses philosopher/MSFragger to process whole cell lysate TMT MS data from CPTAC
# this cancer cell paper used this method (published in 2021, March, 100 citations in 2022, Nov): https://pubmed.ncbi.nlm.nih.gov/33417831/
# https://github.com/Nesvilab/philosopher/wiki/Pipeline-mode-for-TMT-analysis
###################################################
## changes in version 2: (all changes are in the philosopher.yml file)
## 1. have two different philosopher.yml for TMT or itraq. (settings are from FragPipe preset workflow)
## 2. for TMT experiment, change the digest_mass_range from 500-5000 to 200-5000
## 3. for TMTIntegrator, all output using both unique and razor peptide will be generated in the ./unique_and_razor_peptide folder
###################################################
# philosopher workspace --init 
# philosopher database --custom UniProtKB_human_9606_with_Typsin_with_TCGA_candidates_with_HNSCC_candidates.fasta --contam 
## step 1. download data from CPTAC and organized in a specific format
## https://github.com/esacinc/PDC-Public/tree/master/tools/downloadPDCData
## python3 downloadPDCData.py <PDC File manifest.csv>
# python3 /home/yonghao/tricks/downloadPDCData.py  PDC_file_manifest_11032022_114341.csv

## step 2. move mzML files
# ls -d */*/*/* | while read file; do escapespace=${file// /\\\ }; foldername=$(basename "$escapespace") ; echo "mv $escapespace $foldername" >> 1_mv_mzML_files.txt; done
# parallel_GNU -j 20 < 1_mv_mzML_files.txt

## step 3. unzip mzML files
# find . -name *gz | while read file; do echo "gunzip ${file}" >> 2_unzip_mzML_files.txt ; done
# parallel_GNU -j 20 < 2_unzip_mzML_files.txt 

## step 4. generate annotation.txt files in each folder
# python3 /home/yonghao/tricks/CPTAC_metadata_parser.py --input PDC_study_experimental_11072022_105438.csv

## step 5. copy parameter folder
# rsync -avFP /scratch/twlab/yliang/HNSCC/data/CPTAC/test1/params .

## step 6. run philosopher
# how to run this script:
# bash /home/yonghao/tricks/philosopher_TMT_pipe_htcf_v2.sh --input . --study PDC000110_OV_TMT10 --channel 10 --assay tmt --reftag pool
# bash /home/yonghao/tricks/philosopher_TMT_pipe_htcf_v2.sh --input . --study PDC000114_OV_iTRAQ4 --channel 4 --assay itraq --reftag Reference
#########################################################################################################################################################
## INPUT
# Arguments to bash script to decide on various parameters
INPUT_DIR="." # it should have a fastq folder in it
STUDY="PDC000110_OV_TMT10"
CHANNEL_NUM="10"
#ASSAY="tmt" # isobaric labeling brand (tmt, itraq)
REFTAG="pool" # pool for tmt and Reference for itraq (tag to locate the reference channel)

#########################################################################################################################################################

while [[ $# > 1 ]]
do
key="$1"
case $key in
    -i|--input)
    INPUT_DIR="$2"
    shift # past argument
    ;;
    -s|--study)
    STUDY="$2"
    shift # past argument
    ;;
    -c|--channel)
    CHANNEL_NUM="$2"
    shift # past argument
    ;;
    -a|--assay)
    ASSAY="$2"
    shift # past argument
    ;;
    -t|--reftag)
    REFTAG="$2"
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

if [ $ASSAY == "tmt" ]
then
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_tmt.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher.yml
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_uniquepeptide_tmt.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher_uniquepeptide.yml
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_tmt_loose_filter.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher_loose_filter.yml
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_uniquepeptide_tmt_loose_filter.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher_uniquepeptide_loose_filter.yml
elif [ $ASSAY == "itraq" ]
then
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_itraq.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher.yml
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_uniquepeptide_itraq.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher_uniquepeptide.yml
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_itraq_loose_filter.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher_loose_filter.yml
    sed "s|HOLDENCHANNEL|${CHANNEL_NUM}|g" params/philosopher_nochannelnumber_uniquepeptide_itraq_loose_filter.yml | sed "s|HOLDENASSAY|${ASSAY}|g" | sed "s|HOLDENREFTAG|${REFTAG}|g" > params/philosopher_uniquepeptide_loose_filter.yml
else
    echo "please specify tmt or itraq"
fi

## the real issue for TMTIntegrator is that it's reading the yml file from the top directory instead of from the params folder, so I need to put the yml file in the top directory
rsync -avFP params/philosopher.yml .
rsync -avFP params/philosopher_uniquepeptide.yml .
rsync -avFP params/philosopher_loose_filter.yml .
rsync -avFP params/philosopher_uniquepeptide_loose_filter.yml .

mkdir ./unique_and_razor_peptide
mkdir ./unique_and_razor_peptide_loose_filter
mkdir ./unique_peptide_only
mkdir ./unique_peptide_only_loose_filter

DATASETS=`ls | grep -v txt | grep -v params | grep -v csv | grep -v yml | grep -v out | grep -v err | grep -v unique | sed -z 's/\n/ /g'`

jid=`sbatch <<- PHILOSOPHER | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o philosopher_%j_${STUDY}.out
	#SBATCH -e philosopher_%j_${STUDY}.err 
	#SBATCH -c 20
	#SBATCH --mem 135g
	#SBATCH -J "philosopher_${STUDY}"
	date
	echo "bash /home/yonghao/tricks/philosopher_TMT_pipe_htcf_v2.sh --input ${INPUT_DIR} --study ${STUDY} --channel ${CHANNEL_NUM} --assay ${ASSAY} --reftag ${REFTAG}"
	echo "philosopher pipeline --config params/philosopher.yml ${DATASETS}"
	philosopher pipeline --config params/philosopher.yml ${DATASETS}
	date
PHILOSOPHER`
dependphilosopher="afterok:$jid"

## step 7. run TMTIntegrator and only use unique peptide
	jid=`sbatch <<- TMTINTEGRATOR | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o TMTIntegrator-%j_${STUDY}.out
	#SBATCH -e TMTIntegrator-%j_${STUDY}.err 
	#SBATCH -c 5
	#SBATCH --mem 32g
	#SBATCH -J "TMTIntegrator_${STUDY}"
	#SBATCH -d $dependphilosopher    
	date
	echo "java -jar /home/yonghao/softwares/TMTIntegrator_v4.0.2/TMTIntegrator_v4.0.2.jar philosopher_uniquepeptide.yml */psm.tsv"
	java -jar /home/yonghao/softwares/TMTIntegrator_v4.0.2/TMTIntegrator_v4.0.2.jar philosopher_uniquepeptide.yml */psm.tsv
    date
    echo "java -jar /home/yonghao/softwares/TMTIntegrator_v4.0.2/TMTIntegrator_v4.0.2.jar philosopher_loose_filter.yml */psm.tsv"
    java -jar /home/yonghao/softwares/TMTIntegrator_v4.0.2/TMTIntegrator_v4.0.2.jar philosopher_loose_filter.yml */psm.tsv
    date
    echo "java -jar /home/yonghao/softwares/TMTIntegrator_v4.0.2/TMTIntegrator_v4.0.2.jar philosopher_uniquepeptide_loose_filter.yml */psm.tsv"
    java -jar /home/yonghao/softwares/TMTIntegrator_v4.0.2/TMTIntegrator_v4.0.2.jar philosopher_uniquepeptide_loose_filter.yml */psm.tsv
	date
TMTINTEGRATOR`

# rm -rf ./params ./philosopher* ./*xml ./*tsv ./*err ./*out ./unique* */*tsv */interact* */*pepXML */*fas */.meta

