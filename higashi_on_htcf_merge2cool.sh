#!/bin/bash
#!/usr/bin/awk
#####################################################################
# Run higashi on htcf
# /scratch/twlab/yliang/Oocyte_scHi-C/analysis/higashi/flyamer_test1_test2
# interactive -p gpu
# conda activate higashi
# bash /home/yonghao/tricks/higashi_on_htcf.sh --input /scratch/twlab/yliang/Oocyte_scHi-C/analysis/higashi/flyamer_test1_test2 --configuration flyamer_test1_test2_htcf.json
#####################################################################
INPUT_PATH="."
CONFIGURATION="./flyamer_test1_test2_htcf.json"
debugdir="./debug"
OUTPUT="flyamer_SN"
LIST="SN_cell_list.txt"
TYPE="selected"

while [[ $# > 1 ]]
do
	key="$1"
	case $key in
		--input)
		INPUT_PATH="$2"
		shift # pass argument
		;;
		--configuration)
		CONFIGURATION="$2"
		shift # pass argument
		;;
		--debug)
		debugdir="$2"
		shift # pass argument
		;;
		--output)
		OUTPUT="$2"
		shift # pass argument
		;;
		--list)
		LIST="$2"
		shift # pass argument
		;;
		--type)
		TYPE="$2"
		shift # pass argument
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

cd $INPUT_PATH

	jid=`sbatch <<- MERGE | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/merge2cool-%j.out
	#SBATCH -e $debugdir/merge2cool-%j.err 
	#SBATCH -J "3_merge2cool"
	#SBATCH --mem=128G

	date
	python3 /home/yonghao/softwares/higashi/Higashi/higashi/Merge2Cool.py -c $CONFIGURATION -o $OUTPUT -l $LIST -t $TYPE
	date
MERGE`

