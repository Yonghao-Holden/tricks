#!/bin/bash
# programmer: Holden Liang
# this script uses fragpipe/MSFragger to process whole cell lysate TMT MS data from CPTAC
# this cancer cell paper used philosopher (published in 2021, March, 100 citations in 2022, Nov): https://pubmed.ncbi.nlm.nih.gov/33417831/
# which is coming from the same lab and uses MSFragger
# however, philosopher only support regular TMT analysis, but doesn't support open search. And Fragpipe contains pre-built workflow that's easy to use
# SO, instead of using philosopher, i'm switching to Fragpipe
# https://github.com/Nesvilab/FragPipe
###################################################
# philosopher workspace --init 
# philosopher database --custom uniprot_human.fasta --contam 

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

## step 5. generate manifest file for fragpipe run
# python3 /home/yonghao/tricks/fragpipe_manifest_generater.py --study PDC000110_OV_TMT10

## step 6. run fragpipe
# how to run this script:
# bash /home/yonghao/tricks/fragpipe_htcf_v1.sh --input . --study PDC000110_OV_TMT10 --database /scratch/twlab/yliang/HNSCC/data/CPTAC/test1_database/2022-12-05-decoys-contam-uniprot_human.fasta.fas --workflow TMT10 --reftag pool
#########################################################################################################################################################
## INPUT
# Arguments to bash script to decide on various parameters
INPUT_DIR="." 
STUDY="PDC000110_OV_TMT10"
DATABASE="/scratch/twlab/yliang/HNSCC/data/CPTAC/test1_database/2022-12-05-decoys-contam-uniprot_human.fasta.fas"
WORKFLOW="TMT10" # check available workflow in /home/yonghao/softwares/FragPipe_v19.0/fragpipe/workflows
MANIFEST="" # path to manifest file
OUTPUT_DIR="./fragpipe_result"
MEM=128
THREAD=20
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
    --database)
    DATABASE="$2"
    shift # past argument
    ;;
    -w|--workflow)
    WORKFLOW="$2"
    shift # past argument
    ;;
    -m|--manifest)
    MANIFEST="$2"
    shift # past argument
    ;;
    -o|--output)
    OUTPUT_DIR="$2"
    shift # past argument
    ;;
    --mem)
    MEM="$2"
    shift # past argument
    ;;
    --thread)
    THREAD="$2"
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

cd ${INPUT_DIR}

mkdir ${OUTPUT_DIR}

# https://github.com/Nesvilab/FragPipe/blob/3e0189ecb8478ed8ce614bfd8bbcdee98ea05b5c/MSFragger-GUI/src/com/dmtavt/fragpipe/ui.translations
sed "s|tmtintegrator.ref_tag=Bridge|tmtintegrator.ref_tag=${REFTAG}|g" /home/yonghao/softwares/FragPipe_v19.0/fragpipe/workflows/${WORKFLOW}.workflow | sed "s|tmtintegrator.add_Ref=1|tmtintegrator.add_Ref=-1|g" | sed "s|tmtintegrator.groupby=0|tmtintegrator.groupby=-1|g" | sed "s|tmtintegrator.unique_pep=false|tmtintegrator.unique_pep=true|g" > ./${WORKFLOW}.workflow

echo "database.db-path=${DATABASE}" >> ./${WORKFLOW}.workflow

jid=`sbatch <<- PHILOSOPHER | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o fragpipe_%j_${STUDY}.out
	#SBATCH -e fragpipe_%j_${STUDY}.err 
	#SBATCH -c 20
	#SBATCH --mem 135g
	#SBATCH -J "fragpipe_${STUDY}"
	date
	echo "fragpipe --headless --workflow ./${WORKFLOW}.workflow --manifest ${STUDY}_manifest.fp-manifest --workdir ${OUTPUT_DIR} --ram ${MEM} --threads ${THREAD} --config-msfragger /home/yonghao/softwares/MSFragger_3.6/MSFragger-3.6/MSFragger-3.6.jar --config-ionquant /home/yonghao/softwares/IonQuant_1.8.9/IonQuant-1.8.9/IonQuant-1.8.9.jar --config-philosopher /home/yonghao/softwares/Philosopher_4.6.0/philosopher --config-python /home/yonghao/softwares/anaconda3/bin/python"
	fragpipe --headless --workflow ./${WORKFLOW}.workflow --manifest ${STUDY}_manifest.fp-manifest --workdir ${OUTPUT_DIR} --ram ${MEM} --threads ${THREAD} --config-msfragger /home/yonghao/softwares/MSFragger_3.6/MSFragger-3.6/MSFragger-3.6.jar --config-ionquant /home/yonghao/softwares/IonQuant_1.8.9/IonQuant-1.8.9/IonQuant-1.8.9.jar --config-philosopher /home/yonghao/softwares/Philosopher_4.6.0/philosopher --config-python /home/yonghao/softwares/anaconda3/bin/python
	date
PHILOSOPHER`

