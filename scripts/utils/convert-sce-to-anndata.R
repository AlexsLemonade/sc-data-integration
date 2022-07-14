# Script used to convert SCE objects to AnnData objects as HDF5 files

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column
# --sce_dir: Path to folder where SCE objects to be converted are stored, 
#   each file should contain the library ID in the filename and be stored as an RDS file.
#   Typically this is the output from running scpca-downstream-analyses 
# --anndata_output_dir: Path to folder where all AnnData files will be saved as HDF5 files 

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
    opt_str = c("--sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "scpca-downstream-analyses"),
    help = "Path to folder where SCE objects to be converted are stored, 
    each file should contain the library ID in the filename and be stored as an RDS file.
    Typically this is the output from running scpca-downstream-analyses"
  ),
  make_option(
    opt_str = c("--anndata_output_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "anndata", "individual_anndata_objects"),
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
library_id <- library_metadata_df %>%
  dplyr::pull(library_biomaterial_id)

# setup output directory 
if(!dir.exists(opt$anndata_output_dir)){
  dir.create(opt$anndata_output_dir, recursive = TRUE)
}

# Identify SCE files -----------------------------------------------------------

# find SCE files that match library ID 
library_search <- paste(library_id, collapse = "|")
all_library_files <- list.files(opt$sce_dir,
                        pattern = library_search, 
                        recursive = TRUE,
                        full.names = TRUE)
# just include RDS files, otherwise HTML files will also be found 
sce_files <- all_library_files[grep(pattern = ".rds", all_library_files)]

# if the number of sce files is different then the library ID's find the missing files
if(length(sce_files) < length(library_id)){
  
  libraries_found <- stringr::str_extract(sce_files, library_search)
  missing_libraries <- setdiff(library_id, libraries_found)
  
  stop(
    glue::glue(
      "Missing SCE object for {missing_libraries}.
      Make sure that you have run `01-run-downstream-analyses.sh`."
    )
  )
}

# Write H5 ---------------------------------------------------------------------

# create paths to anndata h5 files 
anndata_files <- file.path(opt$anndata_output_dir,
                           paste0(library_id, 
                                  "_anndata.h5"))

# small function to read in sce and export as anndata 
export_anndata <- function(sce_file,
                           anndata_file){
  
  sce <- readr::read_rds(sce_file)
  scpcaTools::sce_to_anndata(sce,anndata_file = anndata_file)
}

# apply export_anndata function to all sce files
purrr::walk2(sce_files, anndata_files, export_anndata)
