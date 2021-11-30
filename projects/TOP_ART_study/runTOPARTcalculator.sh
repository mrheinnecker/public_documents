#!/bin/bash

########################
## VARIABLE DEFINITON ##
########################

#### PART 1: file-selection
## directory containing scripts for PART 1
script_dir_one="$(pwd -P)/PART1_file_selection"
## location of R file containing functions (f)
fncts_one="${script_dir_one}/PART1_fncts.R"

## cohort to be analyzed
COHORT="MASTER"
## date of datafreeze
DATA_FREEZE="20211115144100"

## main project directory (d)
project_dir="/icgc/dkfzlsdf/analysis/hipo/hipo_021"
## path to TOP-ART specific directory
topart_dir="${project_dir}/cohort_analysis/TOP-ART"
## path to sub-directory containing necesarry inputs
input_dir="${topart_dir}/input_and_references"
## path to sub-directory to store intermediate data
intermediate_dir="${topart_dir}/intermediate_data/${COHORT}"

## table with final solutions for ploidy purity combinations and other meta data (r)
PLOIPURITABLE="/icgc/dkfzlsdf/analysis/hipo/hipo_021/results_per_pid_masterReprocessing/metadata_curated/metaMASTER_210607.csv"
## onkostar export table (n)
ONKOSTAR_EXP="${input_dir}/OnkoStarExport/20210715_BRCAness_Abfrage_removeUnwantedLineBreaks.csv"
## assignment of reMASTERed baskets to MASTER cohort (b)
BASKET_ASSI="${input_dir}/OnkoStarExport/reMASTERed_baskets/mapping_ICD_reMASTERed_baskets.csv"

## list of pids for which the input files are collected (p)
pids_list="${input_dir}/PID_list/pids_${COHORT}_${DATA_FREEZE}.txt"
## list of collected files per pid-sample-seqmethod (l)
pidsamseq_list="${intermediate_dir}/pid-sample-seqmethod_results_${DATA_FREEZE}.csv"
## table containing the selected files per pid-sample-seqmethod... main input (t)
selected_files="${intermediate_dir}/selected_files_PIDSAM_${DATA_FREEZE}.tsv"


#### PART 2: pid-wise evaluation
## directory containing scripts for PART 2
script_dir_two="$(pwd -P)/PART2_pid_wise_evaluation"
## location of R file containing functions (f)
fncts_two="${script_dir_two}/PART2_fncts.R"
## location of reference genome with annotated promoter regions of all hg19 genes
hg19_annotation="${input_dir}/gene_lists/promoter_regions_2000bp.tsv"

## directory containing results per PID-sample-seq_method combination
results_per_pid_dir="${topart_dir}/results_per_pid"


#### PART 3: data-integration
## directory containing scripts for PART 3
script_dir_three="$(pwd -P)/PART3_data_integration"
## location of R file containing functions (f)
fncts_three="${script_dir_three}/PART3_fncts.R"

## path to sub-directory containing cohord-wide output data
results_dir="${topart_dir}/results_cohort/${COHORT}"
## table with germline varinants created by the colleagues from Dresden
GERM_INFO="${input_dir}/germlineInfoDresden/2020-10-22-NCT-MASTER_Patients_Variants_Genes_Artefacts_curated.xlsx"
## result table from HRDetect pipeline
HRDETECT_TABLE="${project_dir}/HRDetect/data/outputfiles/hrdetect_cohort_table/hrdetect_output_table_1ceac6e80ed8604fb15593210a2abd5f.tsv"
## TOP-ART scores as annotated during the molecular tumorboard
MTB_INFO="${project_dir}/TOP-ART-Scores.txt"
## list of genes for the germline and somatic criteria of the TOP-ART score
GERMLINE_GENES="${input_dir}/gene_lists/genes_germline_criterion.tsv"
SOMATIC_GENES="${input_dir}/gene_lists/genes_somatic_criterion.tsv"
## subtable of selected files contating those that need to be run because of too old versions
to_run="${intermediate_dir}/selected_samples_to_run_${DATA_FREEZE}.tsv"
## result table containing all data layers
integrated_data="${results_dir}/combined_criteria/integrated_data_${DATA_FREEZE}.tsv"


##########
## MAIN ##
##########

## DEPENDENCIES ##

if [ "$(conda env list | grep topart)" != "" ]; then
  # conda create -n topart -c defaults -c bioconda -c conda-forge "r-base>=4.0.0" r-getopt r-tidyverse r-vcfr
  conda activate topart
else
  module load R/4.1.0
fi


## PROCESSING STEPS ##

#### PART 1: file-selection
## create list of result files per pid-sample-seqmethod from list of pids
use_cores=10 ## max 10
bsub \
  -M "1GB" -R "rusage[mem=1GB]" \
  -n "${use_cores}" -R "span[ptile=${use_cores}]" \
  -W "03:30" \
  -env "all" \
  -J "get_pid-sample-seqmethod_results_${DATA_FREEZE}" \
  -o "${intermediate_dir}/pid-sample-seqmethod_results_${DATA_FREEZE}.err" \
  -K \
bash \
  "${script_dir_one}/pid-sample-seqmethod_results/get_pid-sample-seqmethod_results.sh" \
  "${pids_list}" \
  "/icgc/dkfzlsdf/analysis" "hipo/hipo_021" \
  "${PLOIPURITABLE}" \
  "${pidsamseq_list}" \
  "${use_cores}" \
  "ROWS" "INFO"
## change permission
chmod 775 "${pidsamseq_list}"
chmod 775 "${intermediate_dir}/pid-sample-seqmethod_results_${DATA_FREEZE}.err"

## filter list of result files per pid-sample-seqmethod, convert to table and combine with onkostar export
bsub \
  -M "2GB" -R "rusage[mem=2GB]" \
  -W "0:30" \
  -env "all" \
  -J "TOPART_SelectFiles_${DATA_FREEZE}" \
  -o "${intermediate_dir}/selected_files_PIDSAM_${DATA_FREEZE}.err" \
  -K \
Rscript \
  "${script_dir_one}/pid-sample-seqmethod_results/select_files.R" \
  "-f ${fncts_one}" \
  "-r ${PLOIPURITABLE}" \
  "-b ${BASKET_ASSI}" \
  "-n ${ONKOSTAR_EXP}" \
  "-l ${pidsamseq_list}" \
  "-t ${selected_files}"
## change permission
chmod 775 "${selected_files}"
chmod 775 "${intermediate_dir}/selected_files_PIDSAM_${DATA_FREEZE}.err"


#### PART 2: pid-wise evaluation
Rscript ${script_dir_two}/select_samples_to_run.R -s ${selected_files} -o ${to_run} -v "1.0.0"


#selected_files="/home/m168r/home_extension/test_full_cohort_analysis/intermediate_data/MASTER/selected_files_longer_than_2hrs.tsv"
cat $to_run | while  read pid sample seqm yap cnv hrd car germs germi somsnv somindel bam rna comp notneeded meth rest
do

if [[ $comp == "complete" ]]; then

pid_dir="${results_per_pid_dir}/${seqm}/${pid}"

if [ ! -d "$pid_dir" ]; then
mkdir $pid_dir
fi

sample_dir="${pid_dir}/${sample}"

if [ -d "$sample_dir" ]; then
rm -r $sample_dir
fi

mkdir $sample_dir
pidlog="${sample_dir}/run_logfile.txt"
bsub -W 5:00 -M 10000 -R "rusage[mem=10000]" -J "TOP-ART analysis ${pid}" \
-o $pidlog -env "all" \
"Rscript "${script_dir_two}/top_art_calculator.R" -d $GERMLINE_GENES -e $SOMATIC_GENES -g $hg19_annotation -y $yap -c $cnv -h $hrd -a $car -s $somsnv -i $somindel -b $bam -r $rna -f $fncts_two -m $sample_dir -t $sample -n $meth -x "TRUE" -o "TRUE""
fi

done

#### PART 3: data-integration

bsub -W 0:50 -M 20000 -R "rusage[mem=20000]" \
-J "integration of data tables of TOP-ART cohort" \
-o "${results_dir}/integrated_data_${DATA_FREEZE}.err" -env "all" \
"Rscript "${script_dir_three}/integrate_data.R" -d ${GERMLINE_GENES} -s ${SOMATIC_GENES} -t ${results_per_pid_dir} -i ${selected_files} -g ${GERM_INFO} -m ${MTB_INFO} -h ${HRDETECT_TABLE} -o ${integrated_data} -f ${fncts_three} -v ${DATA_FREEZE} -a ${results_dir}"


chmod 775 "${integrated_data}"
chmod 775 "${variant_table}"
