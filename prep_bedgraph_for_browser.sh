#!/bin/bash
#!/usr/bin/awk
# programmer: Holden Liang
# Run it to pack files for browser view

WORK_PATH=$PWD

while [[ $# > 1 ]]
do
	key="$1"
	case $key in
		--prefix)
		prefix="$2"
		shift # pass argument
		;;
		--suffix)
		suffix="$2"
		shift # pass argument
		;;
    	*)
    	        # unknown option
    	;;
	esac
	shift # past argument or value
done

find ${WORK_PATH} -maxdepth 1 -name "${prefix}*${suffix}" | while read file ; do xbase=$(basename $file) ; sort -k1,1 -k2,2n ${file} | bgzip > ${xbase}.for_browser.sorted.gz ; done ;

find ${WORK_PATH} -maxdepth 1 -name "${prefix}*.for_browser.sorted.gz" | while read file ; do tabix -p bed ${file} ; done ;