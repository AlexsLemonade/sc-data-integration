# Script used to create the metadata file needed to run scpca-downstream-analyses

# This script reads in a metadata file containing the desired libraries that have 
# SCE objects to be pre-processed using scpca-downstream-analyses. Before running
# this script, SCE objects either must be converted from loom files using 
# `scripts/00-convert-loom.R` or synced from S3. Using the provided library IDs, 
# the SCE files are identified and if any library ID does not have a corresponding 
# SCE file an error is thrown, printing out those with missing SCE files.
# If all library ID's have SCE files, then filtering is performed 
# using `scpcaTools::filter_counts()`. Finally a new metadata file is created 
# in the expected format for scpca-downstream-analyses with the sample_id, library_id, 
# filtering_method (set to miQC), and filepath to the filtered SCE. 

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column
# --unfiltered_sce_dir: Path to the folder where all unfiltered SCE objects 
#   are saved locally 
# --filtered_sce_dir: Path to the folder where all filtered SCE objects should be
#   saved locally
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
    opt_str = c("--unfiltered_sce_dir"),
    type = "character",
    default = file.path(project_root, "data", "human_cell_atlas", "sce"),
    help = "path to folder where all unfiltered sce objects are located"
  ),
  make_option(
    opt_str = c("--filtered_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "filtered_sce"),
    help = "path to folder where all filtered sce objects are to be stored"
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

# read in library metadata and grab unfiltered sce file paths 
library_metadata_df <- readr::read_tsv(opt$library_file)
library_id <- library_metadata_df %>%
  dplyr::pull(library_biomaterial_id)

unfiltered_sce_files <- file.path(opt$unfiltered_sce_dir, 
                                  library_metadata_df$sce_unfiltered_files_folder, 
                                  library_metadata_df$unfiltered_sce_filename)

# check that unfiltered sce files exist 
# for any files that don't exist, report the library ID 
file_check <- file.exists(unfiltered_sce_files)
if(any(!file_check)){
  missing_libraries <- paste(library_id[which(!file_check)], 
                             collapse = ",")
  stop(
    glue::glue(
      "Missing unfiltered SCE file for {missing_libraries}.
      Make sure that you have synced the S3 file from S3 or run `scripts/00-convert-loom.R` to
      convert the loom file to SCE if necessary."
    )
  )
}

# create output folders for each tissue group/project for sce results 
# save in the same nested folder structure as the unfiltered SCE files
# first grab folders from the metadata and then combine with the path to the filtered_sce_dir
filtered_sce_output_folders <- unique(library_metadata_df$sce_unfiltered_files_folder)
filtered_sce_output_folders <- file.path(opt$filtered_sce_dir, filtered_sce_output_folders)
for (folder in 1:length(filtered_sce_output_folders)){
  if(!dir.exists(filtered_sce_output_folders[folder])){
    dir.create(filtered_sce_output_folders[folder], recursive = TRUE)
  }
}

# Function to filter SCE objects -----------------------------------------------

#' Read in unfiltered SCE, perform filtering of empty droplets, and save filtered SCE
#' Filtering is performed using `DropletUtils::emptyDropsCellRanger`
#'
#' @param unfiltered_sce_file Path to unfiltered SCE object containing empty droplets
#'   saved as an RDS file. 
#' @param filtered_sce_file Path to save filtered SCE object as an RDS file.
#'
#' @return Filtered SCE object 
read_and_filter_sce <- function(unfiltered_sce_file,
                                filtered_sce_file){
  
  unfiltered_sce <- readr::read_rds(unfiltered_sce_file)
  filtered_sce <- scpcaTools::filter_counts(unfiltered_sce)
  
  readr::write_rds(filtered_sce, filtered_sce_file)
  
  return(filtered_sce)
}

# Filter objects ---------------------------------------------------------------

# add filename for filtered sce file to library metadata
library_metadata_df <- library_metadata_df %>%
  dplyr::mutate(filtered_sce_filename = paste0(library_biomaterial_id, "_filtered_sce.rds"))

# build paths to filtered sce files
filtered_sce_files <- file.path(opt$filtered_sce_dir,
                                library_metadata_df$tissue_group,
                                library_metadata_df$project_name,
                                library_metadata_df$filtered_sce_filename)

# apply function to read in unfiltered sce, filter sce, and save filtered sce
filterd_sce_list <- purrr::map2(.x = unfiltered_sce_files,
                                .y = filtered_sce_files,
                                read_and_filter_sce)

# Update metadata --------------------------------------------------------------

# get relative sce file paths by identifying all sce files 
# with library ID in the sce file name in the provided sce directory 
# first create a regular expression corresponding to processed library IDs
sce_file_search <- paste(library_metadata_df$filtered_sce_filename, collapse="|")
# paths to sce files relative to the root directory
filtered_sce_paths <- list.files(project_root, 
                        pattern = sce_file_search, 
                        recursive = TRUE)

# create search term to grab library IDs from filepaths
library_search <- paste(library_id, collapse="|")

# construct metadata in the format needed for downstream analyses 
downstream_df <- data.frame(filepath = filtered_sce_paths) %>%
  # extract library ID from SCE file path 
  dplyr::mutate(library_id = stringr::str_extract(filepath, library_search)) %>%
  # join with metadata to obtain sample ID
  dplyr::left_join(library_metadata_df, by = c("library_id" = "library_biomaterial_id")) %>%
  # add column for filtering method
  dplyr::mutate(filtering_method = "miQC") %>%
  # select only the columns needed for input to downstream analyses
  dplyr::select(sample_biomaterial_id, library_id, filtering_method, filepath) %>%
  # downstream analyses requires sample_id as column name 
  dplyr::rename(sample_id = sample_biomaterial_id) %>%
  readr::write_tsv(opt$output_metadata)

# save updated library metadata with addition of filtered sce filename
readr::write_tsv(library_metadata_df, opt$library_file)
