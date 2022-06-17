# Script used to create the metadata file needed to run scpca-downstream-analyses

# This script reads in a metadata file containing the desired libraries that have 
# SCE objects to be pre-processed using scpca-downstream-analyses. Before running
# this script, SCE objects either must be converted from loom files using 
# `scripts/00-convert-loom.R` or synced from S3. Using the provided library IDs, 
# the SCE files are identified and if any library ID does not have a corresponding 
# SCE file an error is thrown, printing out those with missing SCE files.
# If all library ID's have SCE files, then a new metadata file is created in the 
# expected format for scpca-downstream-analyses with the sample_id, library_id, 
# filtering_method (set to miQC), and filepath. 

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column
# --full_metadata_file: The path to the metadata file for all libraries. This file must 
#   contain columns for `library_biomaterial_id` and `sample_biomaterial_id`
# --sce_dir: Path to the folder where all SCE objects are saved locally 
# --output_metadata: Path to write metadata file to be used to run 
#   scpca-downstream-analyses

# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
library(magrittr)
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--library_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "hca-processed-libraries.tsv"),
    help = "path to file listing all libraries that are to be converted"
  ),
  make_option(
    opt_str = c("-m", "--full_metadata_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "hca-library-metadata.tsv"),
    help = "path to library metadata file where each row is a library"
  ),
  make_option(
    opt_str = c("--sce_dir"),
    type = "character",
    default = file.path(project_root, "data", "human_cell_atlas" , "sce"),
    help = "path to folder where all output sce objects should be stored"
  ),
  make_option(
    opt_str = c("-o", "--output_metadata"),
    type = "character",
    default = file.path(project_root, "sample-info", "hca-downstream-metadata.tsv"),
    help = "path to write metadata file to be used to run scpca-downstream-analyses"
  )
)


# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# File set up ------------------------------------------------------------------

# checks that provided metadata files exist
if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

if(!file.exists(opt$full_metadata_file)){
  stop("--full_metadata_file provided does not exist.")
}

# grab library IDs that have corresponding SCE files 
library_id <- readr::read_tsv(opt$library_file) %>%
  dplyr::pull(library_biomaterial_id)

# read in full metadata file to use later to grab sample ID
full_metadata_df <- readr::read_tsv(opt$full_metadata_file)

# get relative sce file paths by identifying all sce files 
# with library ID in the sce file name in the provided sce directory 
# first create a regular expression from processed library IDs
library_search <- paste(library_id, collapse="|")
sce_files <- list.files(opt$sce_dir, 
                        pattern = library_search, 
                        recursive = TRUE)

# check that all libraries have SCE files, 
# if not stop and warn that SCE file first must be created or synced from S3
id_check <- stringr::str_detect(sce_files, library_search)

if(!all(id_check)){
  missing_libraries <- paste(library_id[which(!id_check)], collapse = ",")
  stop(
    glue::glue(
      "Missing SCE file for {missing_libraries}. 
      Make sure that you have synced the SCE file from S3 or run `scripts/00-convert-loom.R` to 
      convert the loom file to SCE if necessary.
      "
    )
  )
}

# Create metadata file ---------------------------------------------------------

downstream_df <- data.frame(filepath = sce_files) %>%
  # extract library ID from SCE file path 
  dplyr::mutate(library_id = stringr::str_extract(filepath, library_search)) %>%
  # join with full metadata to obtain sample ID
  dplyr::left_join(full_metadata_df, by = c("library_id" = "library_biomaterial_id")) %>%
  # add column for filtering method
  dplyr::mutate(filtering_method = "miQC") %>%
  # select only the columns needed for input to downstream analyses
  dplyr::select(sample_biomaterial_id, library_id, filtering_method, filepath) %>%
  # downstream analyses requires sample_id as column name 
  dplyr::rename(sample_id = sample_biomaterial_id) %>%
  readr::write_tsv(opt$output_metadata)
