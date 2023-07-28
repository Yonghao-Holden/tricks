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

###################################################
## step 1. download data from CPTAC and organized in a specific format
#############
## I download mzML format as it's smaller and is a more general MS data storage format compared to .raw format which is more specific to different company (platform)
## this paper explains more about formats: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3073315/
#############
## this is how PDC (Proteomic Data Commons) introduction for data downloading. 
## https://proteomic.datacommons.cancer.gov/pdc/faq/Multiple_Files
#############
## compared to the gdc-client utility, this python script will also organize files for us in addition to simply downloading them. 
## so i use this python script to download MS data from PDC instead of using the gdc-client.
## https://github.com/esacinc/PDC-Public/tree/master/tools/downloadPDCData
## python3 downloadPDCData.py <PDC File manifest.csv>
#############
# python3 /home/yonghao/tricks/downloadPDCData.py PDC_file_manifest_PDC000221.csv
# python3 /bar/yliang/tricks/downloadPDCData.py

###################################################
## step 2. move mzML files
#############
# ls -d */*/*/* | while read file; do escapespace=${file// /\\\ }; foldername=$(basename "$escapespace") ; echo "mv $escapespace $foldername" >> 2.1_mv_mzML_files.txt; done
# parallel_GNU -j 20 < 2.1_mv_mzML_files.txt
#############
# remove space from folder's name if needed
# find . -maxdepth 1 -d | while read file; do echo "mv ${file// /\\ } ${file// /_}" >> 2.2_remove_space.txt; done
# parallel_GNU -j 20 < 2.2_remove_space.txt
#############
# find . -name *gz | while read file; do echo "mv ${file/Open Standard/Open\\\ Standard} ${file/Open Standard/}" >> 2.3_mv_mzML_files.txt ; done
# parallel_GNU -j 20 < 2.3_mv_mzML_files.txt
# find . -name *raw | while read file; do echo "mv ${file} ${file/Proprietary/}" >> 2.3_mv_raw_files.txt ; done
# parallel_GNU -j 20 < 2.3_mv_raw_files.txt

###################################################
## step 3. unzip mzML files
# find . -name *gz | while read file; do echo "gunzip ${file}" >> 3_unzip_mzML_files.txt ; done
# parallel_GNU -j 20 < 3_unzip_mzML_files.txt 

## step 4. generate annotation.txt files in each folder
# python3 /home/yonghao/tricks/CPTAC_metadata_parser_v2.py --input S054_CPTAC_HNSCC_TMT11_Label_to_Sample_Mapping_File_JHU_May2020_r1.xlsx --data tmt11

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

if [ ! -d ${OUTPUT_DIR} ]
then
    mkdir ${OUTPUT_DIR}
fi

# https://github.com/Nesvilab/FragPipe/blob/3e0189ecb8478ed8ce614bfd8bbcdee98ea05b5c/MSFragger-GUI/src/com/dmtavt/fragpipe/ui.translations
if [[ $WORKFLOW == *"TMT"*  ||  $WORKFLOW == *"iTRAQ"* ]];
then
    sed "s|tmtintegrator.ref_tag=Bridge|tmtintegrator.ref_tag=${REFTAG}|g" /home/yonghao/softwares/FragPipe_v19.1/fragpipe/workflows/${WORKFLOW}.workflow | sed "s|tmtintegrator.add_Ref=1|tmtintegrator.add_Ref=-1|g" | sed "s|tmtintegrator.groupby=0|tmtintegrator.groupby=-1|g" | sed "s|tmtintegrator.unique_pep=false|tmtintegrator.unique_pep=true|g" > ./${WORKFLOW}.workflow
else
    rsync /home/yonghao/softwares/FragPipe_v19.1/fragpipe/workflows/${WORKFLOW}.workflow .
fi

echo "database.db-path=${DATABASE}" >> ./${WORKFLOW}.workflow

mv ./${WORKFLOW}.workflow ${OUTPUT_DIR}/${WORKFLOW}.workflow

# srun -J interactive -p interactive --mem=135GB --cpus-per-task=20 --pty /bin/bash

jid=`sbatch <<- PHILOSOPHER | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o fragpipe_%j_${STUDY}_${OUTPUT_DIR}.out
	#SBATCH -e fragpipe_%j_${STUDY}_${OUTPUT_DIR}.err 
	#SBATCH -c 20
	#SBATCH --mem 135g
	#SBATCH -J "fragpipe_${STUDY}_${OUTPUT_DIR}"
	date
    export PATH=/home/yonghao/softwares/java_19/jdk-19.0.1/bin:\$PATH
	echo "/home/yonghao/softwares/FragPipe_v19.1/fragpipe/bin/fragpipe --headless --workflow ${OUTPUT_DIR}/${WORKFLOW}.workflow --manifest ${MANIFEST} --workdir ${OUTPUT_DIR} --ram ${MEM} --threads ${THREAD} --config-msfragger /home/yonghao/softwares/MSFragger_3.7/MSFragger-3.7/MSFragger-3.7.jar --config-ionquant /home/yonghao/softwares/IonQuant_1.8.10/IonQuant-1.8.10/IonQuant-1.8.10.jar --config-philosopher /home/yonghao/softwares/Philosopher_4.8.1/philosopher --config-python /home/yonghao/softwares/anaconda3/bin/python"
	/home/yonghao/softwares/FragPipe_v19.1/fragpipe/bin/fragpipe --headless --workflow ${OUTPUT_DIR}/${WORKFLOW}.workflow --manifest ${MANIFEST} --workdir ${OUTPUT_DIR} --ram ${MEM} --threads ${THREAD} --config-msfragger /home/yonghao/softwares/MSFragger_3.7/MSFragger-3.7/MSFragger-3.7.jar --config-ionquant /home/yonghao/softwares/IonQuant_1.8.10/IonQuant-1.8.10/IonQuant-1.8.10.jar --config-philosopher /home/yonghao/softwares/Philosopher_4.8.1/philosopher --config-python /home/yonghao/softwares/anaconda3/bin/python
    mv fragpipe_*_${STUDY}_${OUTPUT_DIR}.out ${OUTPUT_DIR}
    mv fragpipe_*_${STUDY}_${OUTPUT_DIR}.err ${OUTPUT_DIR}
    date
PHILOSOPHER`

# rm -rf ./params ./philosopher* ./*xml ./*tsv ./*err ./*out ./unique* */*tsv */interact* */*pepXML */*fas */.meta ./.meta
