#!/bin/bash
#!/usr/bin/awk
# programmer: Holden Liang
# Run it under the folder for eigenvectors
#####################################################################
# This is a pipeline for A/B compartments calling using Matlab.
# Adapted from Dixon, Jesse R., et al. 2015, by Xianglin Zhang.
# Input is .hic files, output is the eigenvector.
#####################################################################

# Arguments to bash script to decide on various parameters
MAXJOBS=6 #Cores per CPU

HIC_PATH="/scratch/yliang/Epi-therapy/data/hi-C/B36_mega_hic" # path to hic files
HIC_FILE_EXTENSION=".hic"
MATRIX_PATH="/scratch/yliang/Epi-therapy/data/hi-C/B36_mega_juicer_50k_KR_observed_matrix" # path to output matrix files

BIN_SIZE=50000 # resolution(50k or 100k)
BIN_SIZE_n=50k
NORMALIZATION=KR

PIPE_PATH="/bar/yliang/softwares/compartment_call_Bing_pipeline_2015"
JUICER_PATH="/bar/yliang/softwares/juicer/CPU/common"
REF_DIR="/bar/yliang/genomes"
PEAK_PATH="/scratch/yliang/Epi-therapy/data/atac-seq/peaknumber_for_compartment"

echo "parallel job limit: ${MAXJOBS}"
echo "hic file location: ${HIC_PATH}"
echo "resolution/bin size: ${BIN_SIZE}"
echo "normalization: ${NORMALIZATION}"

module load java matlab python3

# dump matrix from hic files
find ${HIC_PATH} -maxdepth 1 -name "*${HIC_FILE_EXTENSION}" | while read file ; do xbase=$(basename $file ${HIC_FILE_EXTENSION}) ; echo "bash ${PIPE_PATH}/1_hic2mat.sh --sample $file --basename $xbase --hicpath $HIC_PATH --matrix_path $MATRIX_PATH --bin_size $BIN_SIZE --bin_size_n $BIN_SIZE_n --normalization $NORMALIZATION --REF_DIR $REF_DIR --juicer_path $JUICER_PATH --pipe_path $PIPE_PATH" >> 1_hic2mat_commands.txt ; done ;

echo "(1/3) dumping matrix from hic files"
parallel_GNU -j $MAXJOBS < 1_hic2mat_commands.txt

java -version >> 1_hic2mat_commands.txt

## Compartment call
find ${HIC_PATH} -maxdepth 1 -name "*${HIC_FILE_EXTENSION}" | while read file ; do xbase=$(basename $file ${HIC_FILE_EXTENSION}) ; echo "bash ${PIPE_PATH}/2_mat2compartment.sh --basename ${xbase} --matrix_path ${MATRIX_PATH} --output_path $PWD --bin_size ${BIN_SIZE} --bin_size_n ${BIN_SIZE_n} --normalization ${NORMALIZATION} --pipe_path ${PIPE_PATH}" >> 2_mat2compartment_commands.txt ; done ;

echo "(2/3) calling compartments"
parallel_GNU -j $MAXJOBS < 2_mat2compartment_commands.txt

## Integrate compartment calls with ATAC-seq result to flip the eigenvalue
find ${HIC_PATH} -maxdepth 1 -name "*${HIC_FILE_EXTENSION}" | while read file ; do xbase=$(basename $file ${HIC_FILE_EXTENSION}) ; echo "bash ${PIPE_PATH}/3_compartmentvsatac.sh --basename ${xbase} --bin_size_n ${BIN_SIZE_n} --peak_path ${PEAK_PATH} --pipe_path ${PIPE_PATH} --normalization ${NORMALIZATION} --REF_DIR ${REF_DIR}" >> 3_compartmentvsatac_commands.txt ; done ;

echo "(3/3) integrating with ATAC-seq result"

printf "chr\ttotal bins\ttotal peaks\t%% of '+' bin in the chromosome\t# of '+' bin in the chromosome\tavarage ATAC peaks in one bin\t%% of '-' bin in the chromosome\t# of '-' bin in the chromosome\tavarage ATAC peaks in one bin\t%% of '+' bin in p arm\tavarage ATAC peaks in one bin\t%% of '-' bin in p arm\tavarage ATAC peaks in one bin\tp_plus_peak-p_minus_peak(average)\t%% of '+' bin in q arm\tavarage ATAC peaks in one bin\t%% of '-' bin in q arm\tavarage ATAC peaks in one bin\tq_plus_peak-q_minus_peak(average)\n"\
	>> 3_compartmentvsatac.summary.txt

bash 3_compartmentvsatac_commands.txt






