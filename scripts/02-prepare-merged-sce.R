# Script used to prepare SCE objects to be integrated 
#
# SCE files are grouped by a specified grouping variable that must be present as 
# a column in the `library_file`. SCE objects for each group are merged, highly 
# variable genes for the merged object are calculated, and then merged objects are 
# stored as RDS files in the `--merged_sce_dir`.

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column
# --grouping_var: Column name present in the library metadata file to use for 
#   grouping SCE objects and merging prior to performing HVG selection.
# --sce_dir: Path to folder where SCE objects to be converted are stored, 
#   each file should contain the library ID in the filename and be stored as an RDS file.
#   Typically this is the output from running scpca-downstream-analyses 
# --merged_sce_dir: Path to folder where all merged SCE objects will be stored 
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

# source helper functions 
source(file.path(project_root, "scripts", "utils", "integration-helpers.R"))

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
    help = "Column name present in the library metadata file to use for grouping SCE objects
    and merging prior to performing HVG selection."
  ),
  make_option(
    opt_str = c("--sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "scpca-downstream-analyses"),
    help = "Path to folder where SCE objects to be converted are stored, 
    each file should contain the library ID in the filename and be stored as an RDS file.
    Typically this is the output from running scpca-downstream-analyses"
  ),
  make_option(
    opt_str = c("--merged_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "merged-sce-objects"),
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

# read in library metadata and grab unfiltered sce file paths 
library_metadata_df <- readr::read_tsv(opt$library_file)
library_ids <- library_metadata_df %>%
  dplyr::pull(library_biomaterial_id)

# check that grouping variable is present
if(!opt$grouping_var %in% colnames(library_metadata_df)){
  stop("Must provide a grouping_var that is a column in the library metadata file.")
}

# setup output directory 
if(!dir.exists(opt$merged_sce_dir)){
  dir.create(opt$merged_sce_dir, recursive = TRUE)
}

# Identify SCE files -----------------------------------------------------------

# find SCE files that match library ID 
library_search <- paste(library_ids, collapse = "|")
all_library_files <- list.files(opt$sce_dir,
                        pattern = library_search, 
                        recursive = TRUE,
                        full.names = TRUE)
# just include RDS files, otherwise HTML files will also be found 
sce_files <- all_library_files[grep(pattern = ".rds", all_library_files, ignore.case = TRUE)]

# if the number of sce files is different then the library ID's find the missing files
if(length(sce_files) < length(library_ids)){
  
  libraries_found <- stringr::str_extract(sce_files, library_search)
  missing_libraries <- setdiff(library_ids, libraries_found)
  
  stop(
    glue::glue(
      "\nMissing SCE object for {missing_libraries}.
      Make sure that you have run `01-run-downstream-analyses.sh`."
    )
  )
}

# Merge by group ---------------------------------------------------------------

# get the library IDs from the SCE file names so that we can name the SCEs in the correct order 
sce_file_names <- stringr::str_extract(sce_files, pattern = library_search)
sce_list <- sce_files %>%
  purrr::set_names(sce_file_names) %>%
  purrr::map(readr::read_rds)

# get the names for each group
group_names <- unlist(unique(library_metadata_df[, opt$grouping_var]))

# create split sce list using the group names 
split_sce_list <- split(sce_list, group_names)

# create a list of merged SCE objects by group
merged_sce_list <- split_sce_list %>%
  purrr::map(combine_sce_objects)

# Subset to HVG ----------------------------------------------------------------


#' Identify variable genes for a merged object and add to metadata
#'
#' @param merged_sce SCE object that has been merged using combine_sce_objects
#'
#' @return merged SCE object with variable genes added to metadata
add_var_genes <- function(merged_sce){
 
  # grab variable genes
  var_genes <- perform_hvg_selection(merged_sce)
  
  # add variable genes to metadata
  metadata(merged_sce)$variable_genes <- var_genes
  
  return(merged_sce)
}


# apply HVG calculation to list of merged SCEs
merged_sce_list <- merged_sce_list %>%
  purrr::map(add_var_genes)

# add PCA and UMAP 

merged_sce_list <- merged_sce_list %>%
  purrr::map( ~ perform_dim_reduction(.x, 
                                      var_genes = metadata(.x)$variable_genes,
                                      pca_type = "multi"))

# Write H5 ---------------------------------------------------------------------

# create paths to merged SCE files
# named with the name of the sce list which corresponds to the grouping variable, not library ID
merged_sce_files <- file.path(opt$merged_sce_dir,
                              paste0(names(merged_sce_list), 
                                     "_merged_sce.rds"))

# export files 
purrr::walk2(merged_sce_list, merged_sce_files, readr::write_rds)
