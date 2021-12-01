#!/bin/bash
# working directory (where TOP-ART scripts are loacted)
wrk_dir="/path/to/working_dir/"
# vcf file containing germline variants
germline_variant_file="/path/to/file.vcf"           
# bam file of dna sequencing
bam_dna_file="/path/to/file.bam"                    
# vcf file containing somatic copynumber variants
somatic_cnv_variant_file="/path/to/file.vcf"        
# list of TOP-ART relevant germline genes
top_art_relevant_germline_genes="/path/to/file.tsv" 
# list of TOP-ART relevant somatic genes 
top_art_relevant_somatic_genes="/path/to/file.tsv"  
# gene annotation of required reference genome
hg19_gene_annotation="/path/to/file.tsv"            
# file conatining relevant information of HRD score
hrd_count_file="/path/to/file.tsv"                  
# vcf file containing somatic indels
somatic_indel_variant_file="/path/to/file.vcf"      
# directory to store output if required (store_output==TRUE)
sample_output_dir="/path/to/dir/"                   
# methylation file if present (optional)
methylation_file="/path/to/file.tsv"                
# if FASLE output will only be printed into console, 
# if TRUE output will be stored in given directory
store_output="TRUE"  
# bam file of rna sequencing
bam_rna_file="/path/to/file.bam" 
# vcf file containing somatic snvs
somatic_snv_variant_file="/path/to/file.vcf" 
# name of used sample
sample_name="patientX_sampleX" 
# if FALSE, only TOP-ART relevant genes will be investigated (much faster), 
# if TRUE, TOP-ARt score will be still be calculated according to TOP-ART genes,
# but the same analysis is performde to all genes given in the gene annotation
analyse_all_genes="TRUE" 
# vcf file containing germline variants
yapsa_mutational_signatures_file="/path/to/file.tsv" 
## load R
module load R/4.1.0
## run TOP-ART calculator on one sample
Rscript "${wrk_dir}/TOPARTcalculator.R" \
  "-a ${germline_variant_file}" \
  "-b ${bam_dna_file}" \
  "-c ${somatic_cnv_variant_file}" \
  "-d ${top_art_relevant_germline_genes}" \
  "-e ${top_art_relevant_somatic_genes}" \
  "-f ${wrk_dir}/TOPARTfncts.R" \
  "-g ${hg19_gene_annotation}" \
  "-h ${hrd_count_file}" \
  "-i ${somatic_indel_variant_file}" \
  "-m ${sample_output_dir}" \
  "-n ${methylation_file}" \
  "-o ${store_output}" \
  "-r ${bam_rna_file}" \
  "-s ${somatic_snv_variant_file}" \
  "-t ${sample_name}" \
  "-x ${analyse_all_genes}" \
  "-y ${yapsa_mutational_signatures_file}"

## samples to run must be in a tsv format... the order chosen here is an 
## example of our TOP-ARt workflow
samples_to_run="/path/to/file.tsv"
## main putput directory... need two sub-folder, one for WGS, and one for WES
main_pidwise_dir="path/to/main/output/"
mkdir "${main_pidwise_dir}/WGS"
mkdir "${main_pidwise_dir}/WES"
# start a rowwise loop (one sample per row of input file)
cat $samples_to_run | while  read pid sample seqm yap cnv hrd car germs germi somsnv somindel bam rna comp meth cutoff
do
  ## only proceed if the necessary files are present (completeness column of input tsv file)
  if [[ $comp == "complete" ]]; then
  ## create directory for PatientID (PID), sublevel of sequencing type (WGS/WES)
    pid_dir="${main_pidwise_dir}/${seqm}/${pid}"
    ## only create directory if it is not yet present (there could alredy be older 
    ## samples of the same patient be evaluated and they should not be overwritten)
    if [ ! -d "$pid_dir" ]; then
      mkdir $pid_dir
    fi
    ## create directory for result of the respective sample inside the PID directory
    sample_dir="${pid_dir}/${sample}"
    ## remove the old one if it is already present (~update with new version)
    if [ -d "$sample_dir" ]; then
      rm -r $sample_dir
    fi
    mkdir $sample_dir
    ## create log_file (required if a HPC is used (recommended))
    log_file="${sample_dir}/run_logfile.txt"
    ## example HPC job submission (5:00h; 10Gb memory request)
    bsub \
      -W 5:00 \
      -M 10000 \
      -R "rusage[mem=10000]" \
      -J "TOP-ART analysis ${pid}" \
      -o $pidlog \
      -env "all" \
      Rscript "${script_dir_two}/top_art_calculator.R" \
        -d $top_art_relevant_germline_genes \
        -e $top_art_relevant_somatic_genes \
        -g $hg19_gene_annotation \
        -y $yap \
        -c $cnv \
        -h $hrd \
        -a $car \
        -s $somsnv \
        -i $somindel \
        -b $bam \
        -r $rna \
        -f "${wrk_dir}/TOPARTfncts.R" \
        -m $sample_dir \
        -t $sample \
        -n $meth \
        -x "TRUE" \
        -o "TRUE"
  fi
done