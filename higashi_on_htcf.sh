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
STAGE="1"
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
		--stage)
		STAGE="$2"
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
mkdir $debugdir

## Process data
	jid=`sbatch <<- PROCESS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p gpu
	#SBATCH --gres=gpu
	#SBATCH -o $debugdir/process-%j.out
	#SBATCH -e $debugdir/process-%j.err 
	#SBATCH -J "1_process_data"
	#SBATCH --mem=128G
	
	date
	nvidia-smi 
	date
	python3 /home/yonghao/softwares/higashi/Higashi/higashi/Process.py -c $CONFIGURATION
	date
PROCESS`

dependprocessing="afterok:$jid"

	jid=`sbatch <<- MAINCELL | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p gpu
	#SBATCH --gres=gpu
	#SBATCH -o $debugdir/maincell-%j.out
	#SBATCH -e $debugdir/maincell-%j.err 
	#SBATCH -J "2_main_cell"
	#SBATCH -d $dependprocessing  
	#SBATCH --mem=128G

	date
	python3 /home/yonghao/softwares/higashi/Higashi/higashi/main_cell.py  -c $CONFIGURATION -s $STAGE
	date
MAINCELL`

