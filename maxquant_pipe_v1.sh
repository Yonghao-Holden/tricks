#!/bin/bash
# programmer: Holden
# bash /bar/yliang/tricks/maxquant_pipe_v1.sh -j 3 --study TSTEA_setting --rawfolder /scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/rawfiles
# bash /bar/yliang/tricks/maxquant_pipe_v1.sh -j 3 --study TSTEA_setting_protein --rawfolder /scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/rawfiles --database /scratch/yliang/HNSCC/analysis/candidates/4_binding_affinity_to_mhc/candidates_proteinseq.fa
# bash /bar/yliang/tricks/maxquant_pipe_v1.sh -j 3 --study TSTEA_setting_mini --rawfolder /scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/rawfiles_mini
# bash /bar/yliang/tricks/maxquant_pipe_v1.sh -j 3 --study CPTAC_setting --rawfolder /scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/rawfiles --database /scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/UniProtKB_human_9606_with_Typsin_and_candidates.fasta --mqpar /scratch/yliang/HNSCC/data/MassSpec_PDC000221_CPTAC-HNSCC-Discovery-Study-Proteome/mqpar_CPTAC_setting_v1.xml
#####################################################################
# This is a pipline to perform initial analysis on raw MS spectrum data using MaxQuant 2.1.3.0
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
echo "MaxQuant version:" >> softwares_version.txt
dotnet /bar/yliang/softwares/maxquant/MaxQuant_2.1.3.0/bin/MaxQuantCmd.exe --version 2>> softwares_version.txt

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

