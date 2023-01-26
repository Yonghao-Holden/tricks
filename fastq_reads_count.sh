#!/bin/bash

#find . -name "*.fastq.gz" | while read file ; do echo "echo "${file}" >> read_count.txt; echo \$(zcat "${file}"|wc -l)/4|bc >> read_count.txt" >> count_commands.txt ; done ;

find . -name "*.fastq.gz" | while read file; do echo "echo -e \"${file}\t\$((\$(zcat ${file}|wc -l)/4))\" >> read_count_temp.txt" >> count_commands.txt ; done ;

#find . -name "*.fastq.gz" | while read file; do read_count=$($(zcat ${file}|wc -l)/4|bc) >> count_commands.txt ; done ;

parallel -j 20 < count_commands.txt

sort read_count_temp.txt -k1,1 > read_count.txt

rm read_count_temp.txt
