#!/bin/bash
# programmer: Holden
#####################################################################
# This is a pipline to perform initial analysis on raw MS spectrum data using MSfragger
# https://github.com/Nesvilab/MSFragger
# https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger
# https://github.com/Nesvilab/philosopher/wiki
# https://github.com/Nesvilab/philosopher/wiki/Simple-Data-Analysis
#####################################################################
# default settings
MAXJOBS=6 #Cores per CPU
RAWFOLDER="/scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/rawfiles_mini/"
INPUTEXTENSION="raw"
STUDY="TSTEA_setting"
DATABASE="/scratch/yliang/HNSCC/analysis/candidates/5_novel_peptides/candidate_peptides.fa"
MQPAR="/scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/mqpar_TSTEA_setting.xml"

# Arguments to bash script to decide on various parameters
while [[ $# > 1 ]]
do
key="$1"
case $key in
	-j|--maxjobs)
	MAXJOBS="$2"
	shift # past argument
	;;
	-r|--rawfolder)
	RAWFOLDER="$2"
	shift # past argument
	;;
	-d|--database)
	DATABASE="$2"
	shift # past argument
	;;
	-s|--study)
	STUDY="$2"
	shift # past argument
	;;
	-m|--mqpar)
	MQPAR="$2"
	shift # past argument
	;;
	-i|--inputextension)
	INPUTEXTENSION="$2"
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

cd ${RAWFOLDER}/..
mkdir ${STUDY}
cd ${STUDY}

echo "parallel job limit: ${MAXJOBS}"
echo "raw files location: ${RAWFOLDER}"
echo "input extension: ${INPUTEXTENSION}"
echo "databse: ${DATABASE}"
echo "study names: ${STUDY}"
echo "maxquant setting: ${MQPAR}"

## save softwares_version information
echo "philosopher version:" >> softwares_version.txt
philosopher version 2>> softwares_version.txt

## step 1. create a workspace
## philosopher workspace --clean
philosopher workspace --init

## step 2. create a database
philosopher database --custom ../database/UniProtKB_human_9606_with_Typsin_and_candidates.fasta --contam 

## step 3. generate MSFragger parameter files
/bar/yliang/softwares/java_19/jdk-19.0.1/bin/java -jar /bar/yliang/softwares/MSFragger_3.5/MSFragger-3.5/MSFragger-3.5.jar --config

## step 4. run MSFragger on data
## can't use softlinks for input files
/bar/yliang/softwares/java_19/jdk-19.0.1/bin/java -Xmx32g -jar /bar/yliang/softwares/MSFragger_3.5/MSFragger-3.5/MSFragger-3.5.jar closed_fragger.params /scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/MSFragger/test/01CPTAC_HNSCC_W_JHU_20190507_LUMOS_f01.raw

mv /lodge/data/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/01CPTAC_HNSCC_W_JHU_20190507_LUMOS_f01.pepXML .
mv /lodge/data/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/01CPTAC_HNSCC_W_JHU_20190507_LUMOS_f01_uncalibrated.mgf .

## step 5. validate peptide hits with PeptideProphet
philosopher peptideprophet --database 2022-11-02-decoys-contam-UniProtKB_human_9606_with_Typsin_and_candidates.fasta.fas --decoy rev_ --ppm --accmass --expectscore --decoyprobs --nonparam 01CPTAC_HNSCC_W_JHU_20190507_LUMOS_f01.pepXML

## step 6. perform protein inference and generate a protXML table
philosopher proteinprophet interact-01CPTAC_HNSCC_W_JHU_20190507_LUMOS_f01.pep.xml 


## step 7. filter and estimate FDR
philosopher filter --sequential --razor --picked --tag rev_ --pepxml interact-01CPTAC_HNSCC_W_JHU_20190507_LUMOS_f01.pep.xml  --protxml interact.prot.xml




philosopher workspace --backup







# MaxQuantCmd.exe --changeFolder="<new mqpar.xml>" "<new folder with fasta files>" "<new folder with raw files>" mqpar.xml
find $RAWFOLDER -maxdepth 1 -name "*${INPUTEXTENSION}" | while read file; do 
    xbase=$(basename $file)
    SHORT_MQPAR=$(basename $MQPAR)
    mkdir ${xbase%.*}; cd ${xbase%.*}
    ln -s $file $xbase 
    sed "s|holden_database|$DATABASE|g" ${MQPAR} | sed "s|holden_rawfile|$PWD/$xbase|g" > ${xbase%.*}_${SHORT_MQPAR}
    cd ..
    echo "cd ${xbase%.*}; date > ${xbase%.*}_maxquant.log; dotnet /bar/yliang/softwares/maxquant/MaxQuant_2.1.3.0/bin/MaxQuantCmd.exe ${xbase%.*}_${SHORT_MQPAR} &>> ${xbase%.*}_maxquant.log 2>> ${xbase%.*}_maxquant.err; date >> ${xbase%.*}_maxquant.log" >> maxquant_command.txt;
done

#/bar/yliang/myapps/bin/parallel_GNU -j ${MAXJOBS} < maxquant_command.txt 

