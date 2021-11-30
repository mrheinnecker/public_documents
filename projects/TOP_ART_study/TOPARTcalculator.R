version <- "1.0.2"


options(readr.num_columns = 0)
options(warn=-1)

suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(getopt))
suppressPackageStartupMessages(require(vcfR))
suppressPackageStartupMessages(require(crayon))
suppressPackageStartupMessages(require(stringi))
suppressPackageStartupMessages(require(feather))

#input <- read_tsv("/icgc/dkfzlsdf/analysis/hipo/hipo_021/cohort_analysis/TOP-ART/intermediate_data/MASTER/selected_files_PIDSAM_20210310160500.tsv") %>%
#select(pid=PID, sample, seq_method, yapsa_file, cnv_file=G1_cnv_file, hrd_file=P2_file,  charger_file=G2_file, germline_snv=G2_snv_file,   germline_indel=G2_indel_file, somatic_snv=G1_snv_file,   somatic_indel=G1_indel_file,   bam_file, rna_bam=rna_file, meth_file)%>%
# filter(pid=="H021-JB7K8N", 
        #sample=="tumor_buffy_coat"
 #       )#_metastasis_buffy_coat_WGS


spec <-
  matrix(
    c("charger_file",       "a", 1, "character",
      "bam_file",           "b", 1, "character",
      "cnv_file",           "c", 1, "character",
      "germline_genes",     "d", 1, "character",
      "somatic_genes",      "e", 1, "character",
      "functions",          "f", 1, "character",
      "ref_gen",            "g", 1, "character",
      "hrd_file",           "h", 1, "character",
      "somatic_indel",      "i", 1, "character",
      "main_output_path",   "m", 1, "character",
      "meth_file",          "n", 1, "character",
      "cohort_analysis",    "o", 1, "character", 
      "rna_bam",            "r", 1, "character", 
      "somatic_snv",        "s", 1, "character",
      "sample",             "t", 1, "character",
      "all_genes",          "x", 1, "character",
      "yapsa_file",         "y", 1, "character"
      ), 
    ncol = 4,
    byrow = TRUE )
opt <- getopt(spec)
opt$version <- version

source(opt$functions)
#source("/home/m168r/projects/top-art-study/TopArtCalculator/fncts.R")
#calculate_topart_score(opt)
main(opt)
