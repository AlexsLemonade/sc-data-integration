# Script used to copy ScPCA data processed using scpca-nf for libraries indicated
# in the `sample-info/scpca-processed-libraries.tsv` file. 
# Specifically, only the `_filtered.rds` files are copied for each library to use 
# as input to the `01-run-downstream-analyses.sh` script and ultimately to use for
# data integration. 
# The library metadata file is used to identify which library files to copy to a local 
# directory for filtered SCE files. All missing files are copied, but if all files 
# are already present, no copying is performed. 

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries whose SCE files should be copied.
#    This file must contain the `library_biomaterial_id`, `sample_biomaterial_id`, and `project_name` column
# --filtered_sce_dir: Path to the folder where all filtered SCE files should be stored locally 
# --scpca_nf_bucket: Bucket on S3 where all outputs from scpca-nf can be found

#load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
suppressPackageStartupMessages({
  library(magrittr)
  library(optparse)
})
source(file.path(project_root, "scripts", "utils", "integration-helpers.R"))


# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--library_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "scpca-processed-libraries.tsv"),
    help = "path to metadata file listing all libraries that are to be converted"
  ),
  make_option(
    opt_str = c("--filtered_sce_dir"),
    type = "character",
    default = file.path(project_root, "data", "scpca", "filtered_sce"),
    help = "path to folder where all sce files should be stored"
  ),
  make_option(
    opt_str = c("--scpca_nf_bucket"),
    type = "character",
    default = "s3://nextflow-ccdl-results/scpca-prod/publish",
    help = "Bucket on s3 where sce data is stored"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# File set up ------------------------------------------------------------------
# checks that provided metadata files exist
if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# create directories if they don't exist 
create_dir(opt$filtered_sce_dir)

# identify datasets to be converted
# want to grab both project and library ids together
library_metadata_df <- readr::read_tsv(opt$library_file)

# construct path to filtered RDS files on S3 to copy to local filtered SCE folder 
s3_filtered_sce_files <- library_metadata_df %>%
  dplyr::mutate(sce_filename = paste0(library_biomaterial_id, "_filtered.rds"),
                filtered_file_path = file.path(scpca_project_id,
                                               sample_biomaterial_id,
                                               sce_filename)) %>%
  dplyr::pull(filtered_file_path)

local_filtered_files <- file.path(opt$filtered_sce_dir,
                                  s3_filtered_sce_files)

# Copy SCE files from S3 -------------------------------------------------------

# copy filtered RDS files to local filtered SCE folder (for any files that don't already exist)
# first create a vector of the files that are present on S3 that don't exist in the local directory yet
missing_sce_files <- s3_filtered_sce_files[which(!file.exists(local_filtered_files))]

# if any are missing, copy only those from AWS
if(length(missing_sce_files) > 0){
  
  aws_includes <- paste("--include '", "*", missing_sce_files, "'", sep = '', collapse = ' ')
  
  # build one sync call to copy all missing sce files 
  sync_call <- paste('aws s3 cp', opt$scpca_nf_bucket, opt$filtered_sce_dir, 
                     '--exclude "*"', aws_includes, '--recursive', sep = " ")
  
  system(sync_call, ignore.stdout = TRUE)
  
} else {
  # if all files are present, print a message that there is nothing to copy.
  message("All filtered files present in the `--library_file` already exist in the provided `--filtered_sce_dir.`
          Nothing to copy.")
}

