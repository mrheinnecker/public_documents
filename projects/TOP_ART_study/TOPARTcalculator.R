## this is TOP-ART calculator main Rscript
version <- "1.0.2"
## important options and frequently used packages are loaded globally
options(readr.num_columns = 0)
options(warn=-1)
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(getopt))
suppressPackageStartupMessages(require(vcfR))
suppressPackageStartupMessages(require(crayon))
suppressPackageStartupMessages(require(stringi))
suppressPackageStartupMessages(require(feather))
## all input parameters taken from the bash console are processed here and 
## stored in the main options list-variable
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
## the required functions are loaded and the main method is started
source(opt$functions)
main(opt)