#!/bin/bash
#!/usr/bin/awk
# programmer: Holden
#####################################################################
# This is a pipline to call loop from .hic files using FAN-C
#####################################################################
# default settings
SERVER="htcf"
MAXJOBS="10"
THREADS="20"
HICFOLDER="."
OUTPUTFOLDER="."
HICEXTENSION=".hic"

RESOLUTION="10" # unit is kb. use 1/5/10 kb resolution, 10kb might be better as it has lower sensitivity, so it won't call loops with low enrichment, and we can remove loops that are far away from the diagnoal line by increasing -o in filtering step. 1kb failed
NORMALIZATION="KR" # this is also the default setting in FAN-C
FINAL_PEAKSIZE="25" # unit is kb. use 10/25kb, not a big difference between these two, pick the default
SETTING_PEAKSIZE="2" # this is calculated as int(FINAL_PEAKSIZE-1bp/RESOLUTION)

MIN_DISTANCE="5" # keep this as 50kb, so loops need to be 50kb away from diagnal. calculated as int(50kb/RESOLUTION)
MIN_OBSERVED="5" # no difference between 5, 10 and 25. Loops got filtered out until we set it as 50, which is too stringent

MERGE_DISTANCE="20000" # decreasing this to 10kb will break one loop into two, which are basically the same loop. so keep 20kb, which is the default setting

# Arguments to bash script to decide on various parameters
while [[ $# > 1 ]]
do
key="$1"
case $key in
    -s|--server)
    SERVER="$2"
    shift # past argument
    ;;
    -j|--maxjobs)
    MAXJOBS="$2"
    shift # past argument
    ;;
    -t|--threads)
    THREADS="$2"
    shift # past argument
    ;;
    -h|--hicfolder)
    HICFOLDER="$2"
    shift # past argument
    ;;
    -o|--outputfolder)
    OUTPUTFOLDER="$2"
    shift # past argument
    ;;
    -e|--hicextension)
    HICEXTENSION="$2"
    shift # past argument
    ;;
    -r|--resolution)
    RESOLUTION="$2"
    shift # past argument
    ;;
    -n|--normalization)
    NORMALIZATION="$2"
    shift # past argument
    ;;
    -fp|--final_peaksize)
    FINAL_PEAKSIZE="$2"
    shift # past argument
    ;;
    -sp|--setting_peaksize)
    SETTING_PEAKSIZE="$2"
    shift # past argument
    ;;
    -mid|--min_distance)
    MIN_DISTANCE="$2"
    shift # past argument
    ;;
    -mo|--min_observed)
    MIN_OBSERVED="$2"
    shift # past argument
    ;;
    -med|--merge_distance)
    MERGE_DISTANCE="$2"
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

cd $OUTPUTFOLDER

echo "parallel job limit: ${MAXJOBS}"
echo ".hic file location: ${HICFOLDER}"
echo "hic extension: ${HICEXTENSION}"

## save softwares_version information
echo "FAN-C version:" >> softwares_version.txt
fanc --version >> softwares_version.txt

mkdir intermediate_files
cd intermediate_files

## step 1. Annotating pixels for loop calling (pixel is 2D bin in the contact map)
if [[ $SERVER == "htcf" ]]
then
    mkdir ../scripts
    for chr in {1..22}; do
        find $HICFOLDER -maxdepth 1 -name "*${HICEXTENSION}" | while read file; do xbase=$(basename $file ${HICEXTENSION}); echo -e '#!/bin/bash\n#SBATCH --mem=56G\n#SBATCH --nodes=1\n#SBATCH --time=10-00:00:00\n' >> ../scripts/${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.sbatch;\
            echo "fanc loops ${file}@${RESOLUTION}kb@${NORMALIZATION} ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops --threads $THREADS --chromosome ${chr} --peak-size ${SETTING_PEAKSIZE} &> ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.log" >> ../scripts/${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.sbatch; done; 
    done
    chr="X"
    find $HICFOLDER -maxdepth 1 -name "*${HICEXTENSION}" | while read file; do xbase=$(basename $file ${HICEXTENSION}); echo -e '#!/bin/bash\n#SBATCH --mem=56G\n#SBATCH --nodes=1\n#SBATCH --time=10-00:00:00\n' >> ../scripts/${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.sbatch;\
        echo "fanc loops ${file}@${RESOLUTION}kb@${NORMALIZATION} ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops --threads $THREADS --chromosome ${chr} --peak-size ${SETTING_PEAKSIZE} &> ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.log" >> ../scripts/${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.sbatch; done; 
    find ../scripts -maxdepth 1 -name "*sbatch" | while read file; do echo "sbatch $file" >> 1_annotate_pixels_commands.txt; done;
    echo "(1/) annotate pixels"
    parallel_GNU -j $MAXJOBS < 1_annotate_pixels_commands.txt
else
    for chr in {1..22}; do
        find $HICFOLDER -maxdepth 1 -name "*${HICEXTENSION}" | while read file ; do xbase=$(basename $file ${HICEXTENSION}); echo "fanc loops ${file}@${RESOLUTION}kb@${NORMALIZATION} ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize{FINAL_PEAKSIZE}k_chr${chr}.loops --threads $THREADS --chromosome ${chr} --peak-size ${SETTING_PEAKSIZE} &> ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.log" >> 1_annotate_pixels_commands.txt ; done ; 
    done
    chr="X"
    find $HICFOLDER -maxdepth 1 -name "*${HICEXTENSION}" | while read file ; do xbase=$(basename $file ${HICEXTENSION}); echo "fanc loops ${file}@${RESOLUTION}kb@${NORMALIZATION} ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize{FINAL_PEAKSIZE}k_chr${chr}.loops --threads $THREADS --chromosome ${chr} --peak-size ${SETTING_PEAKSIZE} &> ${xbase}_${NORMALIZATION}_res${RESOLUTION}kb_peaksize${FINAL_PEAKSIZE}k_chr${chr}.loops.log" >> 1_annotate_pixels_commands.txt ; done ; 
    echo "(1/) annotate pixels"
    parallel_GNU -j 1 < 1_annotate_pixels_commands.txt
fi

if [[ $SERVER != "htcf" ]]
then
    ## step 2. Filtering annotated pixels
    find . -maxdepth 1 -name "*loops" | while read file; do echo "fanc loops $file ${file/.loops/.filtered_d${MIN_DISTANCE}_o${MIN_OBSERVED}.loops} --rh-filter -d ${MIN_DISTANCE} -o ${MIN_OBSERVED} &> ${file/.loops/.filtered_d${MIN_DISTANCE}_o${MIN_OBSERVED}.loops.log}" >> 2_filter_pixels_commands.txt; done;
    echo "(2/) filter pixels"
    parallel_GNU -j $MAXJOBS < 2_filter_pixels_commands.txt
    
    ## step 3. Merging unfiltered pixels into loops
    find . -maxdepth 1 -name "*filtered_d*loops" | while read file; do echo "fanc loops $file ${file/.loops/.merged_md${MERGE_DISTANCE}.loops} -j --remove-singlets --merge-distance ${MERGE_DISTANCE} &> ${file/.loops/.merged_md${MERGE_DISTANCE}.loops.log}" >> 3_merge_pixels_commands.txt; done;
    echo "(3/) merge pixels"
    parallel_GNU -j $MAXJOBS < 3_merge_pixels_commands.txt
    
    ## step 4. Exporting to BEDPE
    find . -maxdepth 1 -name "*merged_md*loops" | while read file; do echo "fanc loops $file -b ${file/.loops/.bedpe} &> ${file/.loops/.bedpe.log}" >> 4_export_to_bedpe_commands.txt; done;
    echo "(4/) export to bedpe"
    parallel_GNU -j $MAXJOBS < 4_export_to_bedpe_commands.txt
    
    ## step 5. merge loops from different chromosome
    cd ..
    for sample in `ls ${HICFOLDER}/*${HICEXTENSION}`; do
        xbase=$(basename $sample ${HICEXTENSION})
        individual_bedpe_files=`ls intermediate_files/$xbase*bedpe`
        cat individual_bedpe_files | awk -v OFS="\t" '{print("chr"$1, $2, $3, "chr"$4, $5, $6, $7, $8, "loop_"NR)}' > $xbase_loops.bedpe
        awk -v OFS="\t" '{print($1,$2,$3,$9"_left",$8); print($4,$5,$6,$9"_right",$8)}' $xbase_loops.bedpe > $xbase_loops.bed
    done
    
    ## step 6. aggregate plot
    find $HICFOLDER -maxdepth 1 -name "*${HICEXTENSION}" | while read file ; do xbase=$(basename $sample ${HICEXTENSION}) ;echo "fanc aggregate --loops --pixels 101 --expected-norm ${file}@${RESOLUTION}kb@${NORMALIZATION} $xbase_loops.bedpe $xbase_loops.loopspixel101expected.agg -p $xbase_loops.loopspixel101expected.pdf" >> 6_aggregate_plot_commands.txt; done
    find $HICFOLDER -maxdepth 1 -name "*${HICEXTENSION}" | while read file ; do xbase=$(basename $sample ${HICEXTENSION}) ;echo "fanc aggregate --loops --pixels 21  --expected-norm ${file}@${RESOLUTION}kb@${NORMALIZATION} $xbase_loops.bedpe $xbase_loops.loopspixel21expected.agg  -p $xbase_loops.loopspixel21expected.pdf" >> 6_aggregate_plot_commands.txt; done
    echo "(6/) aggregate plot"
    parallel_GNU -j $MAXJOBS < 6_aggregate_plot_commands.txt
fi

# test for different parameters
## step 1. Annotating pixels for loop calling (pixel is 2D bin in the contact map)
# default setting for fanc is to use KR normalized matrix for loop calling
# use 10/25 kb peak size --> test 1
# use 5/10 kb resolution  --> test 2 # 10kb might be better as it has lower sensitivity, so it won't call loops with low enrichment, and we can remove loops that are far away from the diagnoal line by increasing -o in filtering step
#fanc loops /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@1kb@KR B36_DMSO_KR_res1kb_peaksize25k.loops --threads 20 --chromosomes 20   &> B36_DMSO_KR_res1kb_peaksize25k.loops.log #--peak-size 24 # it failed
#fanc loops /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@5kb@KR B36_DMSO_KR_res5kb_peaksize25k.loops --threads 30 --chromosomes 20  &> B36_DMSO_KR_res5kb_peaksize25k.loops.log #--peak-size 4 
#fanc loops /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.loops --threads 40 --chromosomes 20  &> B36_DMSO_KR_res10kb_peaksize25k.loops.log # --peak-size 2
#fanc loops /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize10k.loops --threads 20 --chromosomes 20 --peak-size 1 &> B36_DMSO_KR_res10kb_peaksize10k.loops.log
#
### step 2. Filtering annotated pixels
## keep -d as 50kb
## test -o 5, 10, 25, 50 --> test 3
#fanc loops B36_DMSO_KR_res5kb_peaksize25k.loops B36_DMSO_KR_res5kb_peaksize25k.filtered_d10_o5.loops --rh-filter -d 10 -o 5 &> B36_DMSO_KR_res5kb_peaksize25k.filtered_d10_o5.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.loops --rh-filter -d 5 -o 5 &>  B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize10k.loops B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.loops --rh-filter -d 5 -o 5 &>  B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.loops --rh-filter -d 5 -o 10 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.loops --rh-filter -d 5 -o 25 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.loops --rh-filter -d 5 -o 50 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.loops.log
#
### step 3. Merging unfiltered pixels into loops
## when merging, test --merge-distance to 10kb, see if we can get rid of large loop interaction --> test 4
##fanc loops B36_DMSO_KR_res1kb_peaksize25k.filtered_d5_o5.loops B36_DMSO_KR_res1kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops -j --remove-singlets --merge-distance 20000 &> B36_DMSO_KR_res1kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops.log
#fanc loops B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.loops B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops -j --remove-singlets --merge-distance 20000 &> B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.loops  B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.loops -j --remove-singlets &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.loops  B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops -j --remove-singlets --merge-distance 20000 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.loops  B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.loops -j --remove-singlets --merge-distance 20000 &> B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.loops -j --remove-singlets --merge-distance 20000 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.loops -j --remove-singlets --merge-distance 20000 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.loops -j --remove-singlets --merge-distance 20000 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.loops.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.loops  B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.loops -j --remove-singlets --merge-distance 10000 &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.loops.log
#
### step 4. Exporting to BEDPE
##fanc loops B36_DMSO_KR_res1kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops -b B36_DMSO_KR_res1kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe B36_DMSO_KR_res1kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe.log
#fanc loops B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops -b B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe &> B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.loops -b B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.bedpe &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.bedpe.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.loops -b B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe.log
#fanc loops B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.loops -b B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.bedpe &> B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.bedpe.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.loops -b B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.bedpe &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.bedpe.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.loops -b B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.bedpe &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.bedpe.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.loops -b B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.bedpe &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.bedpe.log
#fanc loops B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.loops  -b B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.bedpe &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.bedpe.log
#
### step 6. aggregate plot
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.agg -p B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel100expected.pdf &> B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.agg -p B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel21expected.pdf &> B36_DMSO_KR_res5kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel21expected.pdf.log
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.loopspixel100expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.loopspixel21expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_default.loopspixel21expected.pdf.log
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel100expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel21expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md20kb.loopspixel21expected.pdf.log
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.loopspixel100expected.pdf &> B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.loopspixel21expected.pdf &> B36_DMSO_KR_res10kb_peaksize10k.filtered_d5_o5.merged_md20kb.loopspixel21expected.pdf.log
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.loopspixel100expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.loopspixel21expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o10.merged_md20kb.loopspixel21expected.pdf.log
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.loopspixel100expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.loopspixel21expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o25.merged_md20kb.loopspixel21expected.pdf.log
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.loopspixel100expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.loopspixel21expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o50.merged_md20kb.loopspixel21expected.pdf.log
#fanc aggregate --loops --pixels 100 --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.loopspixel100expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.loopspixel100expected.pdf.log
#fanc aggregate --loops --pixels 21  --expected-norm /scratch/yliang/Epitherapy_3D/data/hic/Mega_inter30/B36_DMSO.hic@10kb@KR B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.bedpe B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.agg -p B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.loopspixel21expected.pdf &> B36_DMSO_KR_res10kb_peaksize25k.filtered_d5_o5.merged_md10kb.loopspixel21expected.pdf.log




