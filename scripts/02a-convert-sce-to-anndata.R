# Script used to convert SCE objects to AnnData objects as HDF5 files

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column
# --grouping_var: Column name present in the library metadata file that was used to 
#   indicate which SCE objects should be integrated in `02-prepare-merged-sce.R`
# --merged_sce_dir: Path to folder where merged SCE objects to be converted are stored, 
#   each file should contain one of the values in the `grouping_var` column in the filename 
#   and be stored as an RDS file.
#   Typically this is the output from running scpca-downstream-analyses followed by 
#   merging SCE objects and identifying HVG. 
# --anndata_output_dir: Path to folder where all AnnData files will be saved as HDF5 files 
# 
# **Note that any columns present in the `rowData` of an SCE object that contains 
# duplicated information, e.g. duplicate gene identifiers, are converted to 
# categorical data by the `anndata` package.
# For context and more information see the comment:
# https://github.com/AlexsLemonade/sc-data-integration/pull/42#issuecomment-1187703050

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
    opt_str = c("-g", "--grouping_var"), 
    type = "character",
    default = "project_name", 
    help = "Column name present in the library metadata file that was used to 
    indicate which SCE objects should be integrated in `02-prepare-merged-sce.R`."
  ),
  make_option(
    opt_str = c("--merged_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "merged-sce-objects"),
    help = "Path to folder where SCE objects to be converted are stored, 
    each file should contain the library ID in the filename and be stored as an RDS file.
    Typically this is the output from running scpca-downstream-analyses"
  ),
  make_option(
    opt_str = c("--anndata_output_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "anndata", "merged_anndata_objects"),
    help = "path to folder where all AnnData files will be saved as HDF5 files"
  )
)

# Setup ------------------------------------------------------------------------

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# checks that provided metadata files exist
if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# read in library metadata and grab names of grouping variable 
library_metadata_df <- readr::read_tsv(opt$library_file)

if(!opt$grouping_var %in% colnames(library_metadata_df)){
  stop("--grouping_var must correspond to a column provided in the library metadata file.")
}

group_names <- library_metadata_df %>%
  dplyr::pull(opt$grouping_var) %>%
  unique()

# setup output directory 
if(!dir.exists(opt$anndata_output_dir)){
  dir.create(opt$anndata_output_dir, recursive = TRUE)
}

# Identify merged SCE files ----------------------------------------------------

# find SCE files that match library ID 
group_name_search <- paste(group_names, collapse = "|")
sce_files <- list.files(opt$merged_sce_dir,
                        pattern = group_name_search, 
                        recursive = TRUE,
                        full.names = TRUE)

# if the number of sce files is different then the number of groups to integrate find the missing files
if(length(sce_files) < length(group_names)){
  
  groups_found <- stringr::str_extract(sce_files, group_name_search)
  missing_groups <- setdiff(group_names, groups_found)
  
  stop(
    glue::glue(
      "Missing SCE object for {missing_groups}.
      Make sure that you have run `02-prepare-merged-sce.R`."
    )
  )
}

# Write H5 ---------------------------------------------------------------------

# extract group names from file path to make sure we are writing the input to the correctly named output 
sce_file_names <- stringr::str_extract(sce_files, pattern = group_name_search)

# create paths to anndata h5 files 
anndata_files <- file.path(opt$anndata_output_dir,
                           paste0(sce_file_names, 
                                  "_anndata.h5"))

# small function to read in sce and export as anndata 
export_anndata <- function(sce_file,
                           anndata_file){
  
  sce <- readr::read_rds(sce_file)
  scpcaTools::sce_to_anndata(sce,anndata_file = anndata_file)
}

# apply export_anndata function to all sce files
purrr::walk2(sce_files, anndata_files, export_anndata)
