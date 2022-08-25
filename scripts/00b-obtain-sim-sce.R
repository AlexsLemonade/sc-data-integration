# Script used to obtain SCE objects from simulated data found in HDF5 files

# This script reads in a metadata file containing the desired libraries to obtain 
# SCE objects. If the RDS files already exist for the desired libraries on S3, 
# they are copied from an S3 bucket to a specified local 
# data directory. For libraries that do not have corresponding SCE objects as 
# RDS files, the hdf5 file will be read in and converted to an SCE object before 
# saving the RDS file. Normalization is also performed if HDF5 files are converted. 
# If the hdf5 file is not present locally, the file will be 
# grabbed from the S3 bucket. All SCE objects will then be synced back to S3. 

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from hdf5 to SCE. This file must contain the 
#   `library_biomaterial_id` and `hdf5_filename` column
# --h5_dir: Path to the folder where all hdf5 files should be stored locally 
# --sce_output_dir: Path to the folder where all SCE objects should be saved locally 
# --s3_h5_bucket: Bucket on S3 where hdf5 data can be found 
# --s3_sce_bucket: Bucket on S3 where SCE objects are stored
# --h5_group_column: Column name in colData of the AnnData/SCE object that will be used to split 
#'   SCE objects. Default is "Batch"
# --copy_s3: indicates whether or not to copy existing SCE file from S3 first. 
#   To copy files use `--copy_s3`
# --overwrite: Indicates whether or not to redo hdf5 to SCE conversion and 
#   overwrite any existing SCE files. To overwrite use `--overwrite`

# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
library(magrittr)
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--library_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "control-processed-libraries.tsv"),
    help = "path to metadata file listing all libraries that are to be converted"
  ),
  make_option(
    opt_str = c("--h5_dir"),
    type = "character",
    default = file.path(project_root, "data", "scib_simulated", "hdf5"),
    help = "path to folder where all loom files should be stored"
  ),
  make_option(
    opt_str = c("--s3_h5_bucket"),
    type = "character",
    default = "s3://sc-data-integration/scib_simulated_data/hdf5",
    help = "Bucket on s3 where loom data is stored"
  ),
  make_option(
    opt_str = c("--sce_output_dir"),
    type = "character",
    default = file.path(project_root, "data", "scib_simulated", "sce"),
    help = "path to folder where all output sce objects should be stored"
  ),
  make_option(
    opt_str = c("--s3_sce_bucket"),
    type = "character",
    default = "s3://sc-data-integration/scib_simulated_data/sce",
    help = "Bucket on s3 where SCE objects are stored"
  ),
  make_option(
    opt_str = c("--h5_group_column"),
    type = "character",
    default = "Batch",
    help = "Column name in colData of the AnnData/SCE object that will be used to split 
     SCE objects. Default is 'Batch'"
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
    help = "indicates whether or not to redo HDF5 to SCE conversion and 
      overwrite any existing SCE files. To overwrite use `--overwrite`"
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
if(!dir.exists(opt$h5_dir)){
  dir.create(opt$loom_dir, recursive = TRUE)
}

if(!dir.exists(opt$sce_output_dir)){
  dir.create(opt$sce_output_dir, recursive = TRUE)
}


# identify datasets to be converted
# will need both library ID and h5 filename
library_metadata_df <- readr::read_tsv(opt$library_file)
library_id <- library_metadata_df %>%
  dplyr::pull(library_biomaterial_id)

# get a metadata with just libraries to be processed
# add file info for sce filepaths
library_metadata_df <- library_metadata_df %>%
  # first make filename for sce file
  dplyr::mutate(local_sce_file = paste0(library_biomaterial_id, "_sce.rds"),
                local_sce_path = file.path(opt$sce_output_dir, local_sce_file))


# Functions for converting hdf5 files -------------------------------------------


#' Subset SCE list, normalize 
#'
#' @param name Name corresponding to object to be subset and normalized from sce_list
#' @param sce_list List of named SCE objects
#' @param sce_file_list List of SCE file paths to use for saving the individual SCE objects
#'   found in the SCE list. The name of the SCE object used in the list must be a part of the 
#'   filename, (e.g. if the library id or part of the library id is in the name of the list, 
#'   the filename can be "library_id_sce.rds")
#'
#' @return Individual SingleCellExperiment object containing logcounts assay
#'
normalize_and_export_sce <- function(name,
                                     sce_list,
                                     sce_file_list){
  
  # extract specific sce from list
  sce <- sce_list[[name]]
  
  
  # find index of filelist that corresponds to name of sce (name of sce file will have library ID)
  sce_file_idx <- grep(name, sce_file_list)
  
  # normalize data
  qclust <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters = qclust)
  sce <- scater::logNormCounts(sce)
  
  # make sure group column is labeled with "celltype" lable 
  sce$celltype <- sce$Group
  
  # save as RDS object with SCE 
  readr::write_rds(sce, sce_file_list[sce_file_idx])
  
  return(sce)
  
}

#' Reads in hdf5 files, converts to SCE objects, normalizes and saves as RDS file
#'
#' @param h5_file Path to hdf5 file containing AnnData object
#' @param sce_file_list List of paths containing `.rds` extension to save individual SCE objects
#' @param h5_group_column Column name in colData of the AnnData/SCE object that will be used to split 
#'   SCE objects
#' 
#' @return list of SingleCellExperiment objects
#'
hdf5_to_sce <- function(h5_file,
                        sce_file_list){
  
  
  # read in H5 file as SCE
  # save X assay as counts
  sce <- zellkonverter::readH5AD(h5_file,
                                 X_name = "counts")
  
  # check that group column is present in colData of SCE 
  if (!h5_group_column %in% colnames(colData(sce))){
    stop("Can't split SCE objects in HDF5 file by specified group. `h5_group_column` not found
         in columns of colData.")
  }
  
  # split sce object by batch
  sce_list <- split(sce, sce[[h5_group_column]])
  
  # normalize and export sces, applying the function across the name of SCEs
  purrr::map(names(sce_list),
             ~ normalize_and_export_sce(name = .x,
                                        sce_list,
                                        sce_file_list))
  
  return(sce_list)
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

# grab hdf5 and convert to SCE -------------------------------------------------

# get list of all expected sce paths 
local_sce_paths <- library_metadata_df %>%
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
  
  # first need to see if the hdf5 file is present to create the SCE 
  # if the hdf5 file isn't there, grab it from AWS before converting
  
  # list of all hdf5 files that correspond to missing sce files
  # obtain the path relative to the hdf5 directory
  hdf5_missing_sce <- library_metadata_df %>%
    dplyr::filter(local_sce_path %in% missing_sce_files) %>%
    dplyr::mutate(missing_filepath = file.path(hdf5_filename)) %>%
    dplyr::pull(missing_filepath)
  
  # construct full paths for searching if hdf5 file exists 
  full_missing_hdf5_path <- file.path(opt$h5_dir,
                                      hdf5_missing_sce)
  
  # if any hdf5 files don't exist for missing SCE's then grab those from AWS S3
  if(!all(file.exists(full_missing_hdf5_path))){
    
    # get list of hdf5 files inside s3 directory to include in copying
    aws_local_copy <- hdf5_missing_sce[which(!file.exists(full_missing_hdf5_path))]
    aws_includes <- paste("--include '", "*", aws_local_copy, "'", sep = '', collapse = ' ')
    
    # build one sync call to copy all missing hdf5 files 
    sync_call <- paste('aws s3 cp', opt$s3_h5_bucket, ".", 
                       '--exclude "*"', aws_includes, '--recursive', sep = " ")
    
    system(sync_call, ignore.stdout = TRUE)
  }
  
  # group metadata by hdf5 filename to get list of individual sce files that should be produced for each H5 file
  grouped_metadata_df <- library_metadata_df %>%
    dplyr::filter(hdf5_filename %in% hdf5_missing_sce) %>%
    dplyr::group_by(hdf5_filename) %>%
    dplyr::summarise(sce_files = list(local_sce_path))
  
  # convert to sce objects and write files 
  sce_list <- purrr::map2(file.path(opt$h5_dir, grouped_metadata_df$hdf5_filename),
                          grouped_metadata_df$sce_files, 
                          hdf5_to_sce(h5_group_column = opt$h5_group_column))
  
  # sync sce output to S3 
  all_sce_files <- unique(library_metadata_df$local_sce_file)
  aws_includes <- paste("--include '", all_sce_files, "'", sep = '', collapse = ' ')
  sync_call <- paste('aws s3 sync', opt$sce_output_dir, opt$s3_sce_bucket, 
                     '--exclude "*"', aws_includes, sep = " ")
  system(sync_call, ignore.stdout = TRUE)
  
}
