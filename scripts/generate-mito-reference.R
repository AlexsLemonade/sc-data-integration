# Script used to create mitochondrial gene list from Gencode GTF file 


# import libraries
library(magrittr)
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--s3_ref_bucket"),
    type = "character",
    default = "s3://sc-data-integration/reference-files/gencode_v27",
    help = ""
  ),
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    default = "human.gencode.v27.annotation.gtf",
    help = ""
  ),
  make_option(
    opt_str = c("--local_gtf_dir"),
    type = "character",
    default = "../reference-files/gtf",
    help = ""
  ),
  make_option(
    opt_str = c("--mito_output"),
    type = "character",
    default = "../reference-files/human.gencode.v27.mitogenes.txt",
    help = ""
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check if gtf file has a local copy 
local_gtf_file <- file.path(opt$local_gtf_dir, opt$gtf_file)
if(!file.exists(local_gtf_file)){
  
  # if file doesn't exist, copy from S3 
  includes <- paste("--include '", "*", opt$gtf_file,"'", "*", sep = '', collapse = ' ')
  sync_call <- paste('aws s3 cp', opt$s3_ref_bucket, opt$local_gtf_dir, 
                     '--exclude "*"', includes, '--recursive', sep = " ")
  system(sync_call, ignore.stdout = TRUE)
  
}

# import gtf file and grab mitochondrial genes 
gtf <- rtracklayer::import(gtf_file)
mitogenes <- gtf[seqnames(gtf) == 'chrM']

# save mitochondrial ensemble gene ids
writeLines(unique(mitogenes$gene_id), mito_out)
