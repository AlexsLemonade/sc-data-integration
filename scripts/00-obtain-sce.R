# Script used to obtain SCE objects

# This script reads in a metadata file containing the desired libraries to obtain 
# SCE objects. If the RDS files already exist for the desired libraries on S3, 
# they are copied from an S3 bucket to a specified local 
# data directory. For libraries that do not have corresponding SCE objects as 
# RDS files, the loom file will be read in and converted to an SCE object before 
# saving the RDS file. If the loom file is not present locally, the file will be 
# grabbed from the S3 bucket. All SCE objects will then be synced back to S3. 
# An updated metadata file will be returned with the SCE file information.

# Option descriptions: 
# 
# --metadata_file: The path to the metadata file for all libraries. This file must 
#   contain columns for `library_biomaterial_id`, `tissue_group`, `project_name`,
#   and `loom_file`
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column
# --loom_dir: Path to the folder where all loom files should be stored locally 
# --sce_output_dir: Path to the folder where all SCE objects should be saved locally 
# --s3_loom_bucket: Bucket on S3 where loom data can be found 
# --s3_sce_bucket: Bucket on S3 where SCE objects are stored
# --copy_s3: indicates whether or not to copy existing SCE file from S3 first. 
#   To copy files use `--copy_s3`
# --overwrite: Indicates whether or not to redo loom to SCE conversion and 
#   overwrite any existing SCE files. To overwrite use `--overwrite`

# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
library(magrittr)
library(optparse)
library(LoomExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-m", "--metadata_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "hca-library-metadata.tsv"),
    help = "path to library metadata file where each row is a library"
  ),
  make_option(
    opt_str = c("-l", "--library_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "hca-processed-libraries.tsv"),
    help = "path to file listing all libraries that are to be converted"
  ),
  make_option(
    opt_str = c("--loom_dir"),
    type = "character",
    default = file.path(project_root, "data", "human_cell_atlas", "loom"),
    help = "path to folder where all loom files should be stored"
  ),
  make_option(
    opt_str = c("--sce_output_dir"),
    type = "character",
    default = file.path(project_dir, "data", "human_cell_atlas", "sce"),
    help = "path to folder where all output sce objects should be stored"
  ),
  make_option(
    opt_str = c("--s3_loom_bucket"),
    type = "character",
    default = "s3://sc-data-integration/human_cell_atlas_data/loom",
    help = "Bucket on s3 where loom data is stored"
  ),
  make_option(
    opt_str = c("--s3_sce_bucket"),
    type = "character",
    default = "s3://sc-data-integration/human_cell_atlas_data/sce",
    help = "Bucket on s3 where SCE objects are stored"
  ),
  optparse::make_option(
    c("-c", "--copy_s3"),
    action = "store_true",
    help = "indicates whether or not to copy existing SCE file from S3 first. 
    To copy files use `--copy_s3`"
  ),
  optparse::make_option(
    c("-o", "--overwrite"),
    action = "store_true",
    help = "indicates whether or not to redo loom to SCE conversion and 
      overwrite any existing SCE files. To overwrite use `--overwrite`"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# File set up ------------------------------------------------------------------

# checks that provided metadata files exist
if(!file.exists(opt$metadata_file)){
  stop("--metadata_file provided does not exist.")
}

if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# create directories if they don't exist 
if(!dir.exists(opt$loom_dir)){
  dir.create(opt$loom_dir, recursive = TRUE)
}

if(!dir.exists(opt$sce_output_dir)){
  dir.create(opt$sce_output_dir, recursive = TRUE)
}

# read in metadata file 
all_metadata_df <- readr::read_tsv(opt$metadata_file)

# identify datasets to be converted
# read in library file and find intersection between metadata and library file
# will need both library ID and loom filename
library_metadata_df <- readr::read_tsv(opt$library_file)
library_id <- library_metadata_df %>%
  dplyr::pull(library_biomaterial_id)
loom_files <- library_metadata_df %>%
  dplyr::pull(loom_filename)

# get a metadata with just libraries to be processed
# add file info for sce filepaths
process_metadata_df <- all_metadata_df %>%
  # use filename as unique identifier to select library file
  dplyr::filter(loom_filename %in% loom_files) %>%
  # first make filename for sce file
  dplyr::mutate(local_sce_file = file.path(tissue_group,
                                           project_name,
                                           paste0(library_biomaterial_id, "_sce.rds")),
                # get path to output folder needed for creating directories 
                output_folder = file.path(opt$sce_output_dir,
                                          tissue_group,
                                          project_name),
                # then make complete path to sce file
                local_sce_path = file.path(opt$sce_output_dir,
                                           local_sce_file))

# create output folders for each tissue group/project for sce results 
sce_output_folders <- unique(process_metadata_df$output_folder)
for (folder in 1:length(sce_output_folders)){
  if(!dir.exists(sce_output_folders[folder])){
    dir.create(sce_output_folders[folder], recursive = TRUE)
  }
}


# Function for converting loom files -------------------------------------------

#' Reads in loom files, converts to SCE objects and saves as RDS file
#'
#' @param loom_file Path to loom file containing SingleCellLoomExperiment
#' @param sce_file Path containing `.rds` extension to save SCE object
#'
#' @return SingleCellExperiment object
#'
loom_to_sce <- function(loom_file,
                        sce_file){
  
 
  # import loom file as SingleCellLoomExperiment
  loom <- import(loom_file, type= "SingleCellLoomExperiment")
  
  # first grab matrix from loom file 
  # convert to sparse matrix from DelayedMatrix, sparse is required for SCE
  mtx <- loom@assays@data$matrix
  counts <- as(mtx, "dgCMatrix")
  
  # create new SCE using pieces from loom file 
  sce <- SingleCellExperiment(list(counts = counts),
                              colData = colData(loom),
                              rowData = rowData(loom),
                              metadata = metadata(loom)) 
  
  # add cell barcodes and gene names to column names and row names of SCE
  colnames(sce) = sce$CellID
  rownames(sce) = rowData(sce)$ensembl_ids
  
  # save sce file 
  readr::write_rds(sce, sce_file)
  
  return(sce)
}

# Grab existing SCE from S3 ----------------------------------------------------

# only copy from S3 if option to copy is used at the command line
if(!is.null(opt$copy_s3)){
  # grab all library ID's that should have SCE's copied over
  libraries_include <- paste("--include '", "*", library_id,"'", "*", sep = '', collapse = ' ')
  sync_call <- paste('aws s3 cp', opt$s3_sce_bucket, opt$sce_output_dir, 
                     '--exclude "*"', libraries_include, '--recursive', sep = " ")
  system(sync_call, ignore.stdout = TRUE) 
}

# grab loom and convert to SCE --------------------------------------------------

# get list of all expected sce paths 
local_sce_paths <- process_metadata_df %>%
  dplyr::pull(local_sce_path)

# create a list of sce files that need to be created 
# if overwriting existing SCEs grab all sce paths 
if(!is.null(opt$overwrite)){
  missing_sce_files <- local_sce_paths
} else {
  # if not overwriting existing ones, find which SCE's don't exist yet
  missing_sce_files <- local_sce_paths[which(!file.exists(local_sce_paths))]  
}


# if any libraries are missing a corresponding sce, then create that sce 
if(length(missing_sce_files) != 0){
  
  # first need to see if the loom file is present to create the SCE 
  # if the loom file isn't there, grab it from AWS before converting
  
  # list of all loom files that correspond to missing sce files
  # obtain the path relative to the loom directory
  loom_missing_sce <- process_metadata_df %>%
    dplyr::filter(local_sce_path %in% missing_sce_files) %>%
    dplyr::mutate(missing_filepath = file.path(tissue_group,
                                               project_name,
                                               bundle_uuid,
                                               loom_filename)) %>%
    dplyr::pull(missing_filepath)
  
  # construct full paths for searching if loom file exists 
  full_missing_loom_path <- file.path(opt$loom_dir,
                                      loom_missing_sce)
    
  # if any loom files don't exist for missing SCE's then grab those from AWS S3
  if(!all(file.exists(full_missing_loom_path))){
    
    # get list of loom files inside s3 directory to include in copying
    aws_local_copy <- loom_missing_sce[which(!file.exists(full_missing_loom_path))]
    aws_includes <- paste("--include '", "*", aws_local_copy, "'", sep = '', collapse = ' ')
    
    # build one sync call to copy all missing loom files 
    sync_call <- paste('aws s3 cp', opt$s3_loom_bucket, ".", 
                       '--exclude "*"', aws_includes, '--recursive', sep = " ")
    
    system(sync_call, ignore.stdout = TRUE)
  }
  
  # convert to sce objects and write files 
  sce_list <- purrr::map2(loom_file_paths,
                          missing_sce_files, 
                          loom_to_sce)
  
  # sync sce output to S3 
  all_sce_files <- unique(process_metadata_df$local_sce_file)
  aws_includes <- paste("--include '", all_sce_files, "'", sep = '', collapse = ' ')
  sync_call <- paste('aws s3 sync', opt$sce_output_dir, opt$s3_sce_bucket, 
                     '--exclude "*"', aws_includes, sep = " ")
  system(sync_call, ignore.stdout = TRUE)
  
  # create data frame with sce filenames to add to existing library metadata
  library_search <- paste(library_id, collapse="|")
  sce_filepath_df <- data.frame(unfiltered_sce_filename = missing_sce_files) %>%
    dplyr::mutate(library_biomaterial_id = stringr::str_extract(unfiltered_sce_filename, library_search))
  
  # modify existing library metadata with sce file information
  library_metadata_updated <- library_metadata_df %>%
    dplyr::left_join(sce_filepath_df) %>%
    dplyr::mutate(sce_unfiltered_files_s3_bucket = opt$s3_sce_bucket,
                  sce_unfiltered_files_folder = file.path(tissue_group, project_name),
                  unfiltered_sce_filename = stringr::word(unfiltered_sce_filename, -1, sep = "/")) %>%
    # write out updated library file with new sce file information
    readr::write_tsv(opt$library_file)
  
}

