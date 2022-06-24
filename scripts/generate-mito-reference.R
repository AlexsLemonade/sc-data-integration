# Script used to create mitochondrial gene list from Gencode GTF file 

# First the script searches for the gtf file locally and if not present, 
# will copy the gtf file from the provided S3 bucket to a local directory. 
# The mitochondrial genes are then extracted from the Gencode GTF file and saved 
# to a txt file. 
# Note: This only works for gtf files from Gencode as it identifies mitochondrial 
# genes using `ChrM` in contrast to the chromosome annotation used in Ensembl, `MT`

# Option descriptions: 
#
# --s3_ref_bucket: S3 bucket where gtf file is stored and can be copied from 
# --gtf_file: File name for gtf file 
# --local_gtf_dir: Local directory where gtf file exists or should be copied to.
#   This directory should be ignored by git.
# --mito_output: Path to file to save mitochondrial gene list extracted from 
#   input gtf file.

# import libraries
library(magrittr)
library(optparse)
library(GenomicFeatures)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--s3_ref_bucket"),
    type = "character",
    default = "s3://sc-data-integration/reference-files/gencode_v27",
    help = "S3 bucket where gtf file is stored"
  ),
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    default = "gencode.v27.annotation.gtf",
    help = "File name for gtf file"
  ),
  make_option(
    opt_str = c("--local_gtf_dir"),
    type = "character",
    default = "../reference-files/gtf",
    help = "Local directory where gtf file exists or should be copied to"
  ),
  make_option(
    opt_str = c("--mito_output"),
    type = "character",
    default = "../reference-files/gencode.v27.mitogenes.txt",
    help = "Path to file to save mitochondrial gene list extracted 
      from input gtf file."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Copy GTF File ----------------------------------------------------------------

# check if gtf file has a local copy 
local_gtf_file <- file.path(opt$local_gtf_dir, opt$gtf_file)
if(!file.exists(local_gtf_file)){
  
  # if file doesn't exist, copy from S3 
  includes <- paste("--include '", "*", opt$gtf_file,"'", "*", sep = '', collapse = ' ')
  sync_call <- paste('aws s3 cp', opt$s3_ref_bucket, opt$local_gtf_dir, 
                     '--exclude "*"', includes, '--recursive', sep = " ")
  system(sync_call, ignore.stdout = TRUE)
  
}

# Extract Mito Genes -----------------------------------------------------------

# import gtf file and grab mitochondrial genes 
gtf <- rtracklayer::import(local_gtf_file)
mitogenes <- gtf[seqnames(gtf) == 'chrM']

# save mitochondrial ensemble gene ids
writeLines(unique(mitogenes$gene_id), opt$mito_output)
