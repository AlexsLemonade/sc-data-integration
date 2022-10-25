# Script used to add celltype information to the processed SCE objects from the 
# Dyer RMS dataset (SCPCP000005)

# The processed SCE data should have been obtained from scpca-nf originally and then 
# run through scpca-downstream-analyses using `01-run-downstream-analyses.sh`

# The celltype information is found in the Seurat objects that are stored on S3 in 
# `s3://sc-data-integration/scpca/rms_dyer_seurat`. These first must be copied to the 
# local data directory using `aws s3 sync`

# There are both individual Seurat objects and integrated Seurat objects. The individual 
# Seurat objects contain a column, `cell.cluster.ids` to annotate cells as tumor or non-tumor 
# cell types. The integrated objects contain an additional column, `annot.clusters` that contain
# the tumor cell subtype. This script combines the results in both of those columns in one singular 
# `celltype` column classifying cells as tumor/non-tumor and by tumor subtype

# Note: Not all SCPCP000005 libraries present and submitted to ScPCA contain cell type information
# If a seurat object is missing and therefore there is no cell type information, the `celltype`
# column is populated with `NA`

# Option descriptions: 
# 
# --library_file: The full path to the file listing all libraries that celltype information should be added to.
# --processed_sce_dir: Full path to folder where all filtered and processed sce files are found, output from 
#   running scpca-downstream-analyses 
# --celltype_sce_dir: Full path to folder where all sce files containing celltypes are stored
# --rms_seurat_dir: Full path to folder storing local copy of Seurat objects for RMS data

# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
suppressPackageStartupMessages({
  library(magrittr)
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
})

source(file.path(project_root, "scripts", "utils", "rms-celltype-helper-functions.R"))

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--library_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "rms-processed-libraries.tsv"),
    help = "path to metadata file listing all libraries that celltype information should be added to.
      This file must correspond to the RMS data from the Dyer/Chen project."
  ),
  make_option(
    opt_str = c("--processed_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "scpca", "scpca-downstream-analyses"),
    help = "path to folder where all filtered and processed sce files are found, output from running scpca-downstream-analyses"
  ),
  make_option(
    opt_str = c("--celltype_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "scpca", "celltype_sce"),
    help = "path to folder where all sce files containing celltypes are stored"
  ),
  make_option(
    opt_str = c("--rms_seurat_dir"),
    type = "character",
    default = file.path(project_root, "data", "scpca", "rms_dyer_seurat"),
    help = "path to folder storing local copy of Seurat objects for RMS data"
  )
)


# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# File set up ------------------------------------------------------------------
# checks that provided metadata files exist
if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# check that celltype sce dir exists 
if(!dir.exists(opt$celltype_sce_dir)){
  dir.create(opt$celltype_sce_dir, recursive = TRUE)
}

# check that integrated RMS data is present, needed for grabbing tumor sub-celltypes 
integrated_sn_obj_file <- file.path(opt$rms_seurat_dir, "2020330 RMS integrated snRNA Seurat object-005.Rds")
integrated_sc_pt_obj_file <- file.path(opt$rms_seurat_dir, "20191230 RMS pt tumor non-malignant.Rds")
integrated_pdx_obj_file <- file.path(opt$rms_seurat_dir, "20200818 all PDX malignant cells Seurat object.Rds")

integrated_sce_file_list <- list(
  integrated_sn_obj_file,
  integrated_sc_pt_obj_file,
  integrated_pdx_obj_file
)

# create a combined dataframe with all coldata from integrated objects 
all_integrated_coldata <- integrated_sce_file_list %>%
  purrr::map(grab_integrated_metadata) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(barcode = stringr::word(unique_cell_id, -1, sep = "-")) # add barcode column to coldata

# identify submitter ids and matching sample, library ids
# want to grab both project and library ids together
library_metadata_df <- readr::read_tsv(opt$library_file) %>%
  dplyr::select(sample_biomaterial_id, library_biomaterial_id, submitter_id, project_name, seq_unit) %>%
  dplyr::mutate(sce_processed_filename = paste0(library_biomaterial_id, "_miQC_processed_sce.rds"), 
                sce_processed_filepath = file.path(opt$processed_sce_dir, 
                                                  sample_biomaterial_id, 
                                                  sce_processed_filename),
                seq_unit = dplyr::case_when(seq_unit == "cell" ~ "sc",
                                            seq_unit == "nucleus" ~ "sn"),
                seurat_filename = paste0(submitter_id, "_", seq_unit, ".Rds"),
                seurat_filepath = file.path(opt$rms_seurat_dir,
                                            seurat_filename)) 

# check that all sce filtered objects are present
sce_processed_files <- library_metadata_df %>%
  dplyr::pull(sce_processed_filepath)
missing_sce_files <- sce_processed_files[which(!file.exists(sce_processed_files))]

# if files are missing, grab from scpca-nf bucket 
if(length(missing_sce_files) > 0){
  
  missing_libraries <- basename(missing_sce_files) %>%
    stringr::word(1, sep = "_")
  
  stop(glue::glue("
             Processed SCE files not found for: {missing_libraries}
             
             "))
  
}

# get celltype data for each entry in library metadata from scpca
library_metadata_df %>%
  dplyr::select(sce_processed_filepath, seurat_filepath, library_biomaterial_id, submitter_id) %>%
  purrr::pwalk(add_celltype,
               all_integrated_coldata = all_integrated_coldata)

