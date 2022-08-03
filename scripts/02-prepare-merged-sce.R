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
#   grouping SCE objects and merging.
# --subset_hvg: Indicates whether or not to subset the merged SCE object by highly variable genes.
#   If --subset_hvg is used, the merged SCE object will only contain genes 
#   identified as highly variable genes.
# --use_all_genes: Indicates whether or not to use the all genes as input to performing 
#   principal component analysis. Otherwise only highly variable genes are used 
#   as input.
# --num_genes: Number of highly variable genes to use if --select_hvg is being used.
# --sce_dir: Path to folder where SCE objects to be converted are stored, 
#   each file should contain the library ID in the filename and be stored as an RDS file.
#   Typically this is the output from running scpca-downstream-analyses 
# --merged_sce_dir: Path to folder where all merged SCE objects will be stored 
# 

# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
library(magrittr)
library(optparse)
library(SingleCellExperiment)

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
    and merging."
  ),
  make_option(
    opt_str = c("--subset_hvg"), 
    default = FALSE,
    action = "store_true",
    help = "Indicates whether or not to subset the merged SCE object by highly variable genes.
      If --subset_hvg is used, the merged SCE object will only contain genes 
      identified as highly variable genes."
  ),
  make_option(
    opt_str = c("--use_all_genes"), 
    default = FALSE,
    action = "store_true",
    help = "Indicates whether or not to use the all genes as input to performing 
      principal component analysis. Otherwise only highly variable genes are used 
      as input."
  ),
  make_option(
    opt_str = c("-n", "--num_genes"),
    type = "integer",
    default = 5000,
    help = "Number of highly variable genes to select."
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
    default = file.path(project_root, "results", "human_cell_atlas", "merged_sce"),
    help = "path to folder where all merged SCE objects files will be saved as RDS files"
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

# check that num genes provided is an integer 
if(!is.integer(opt$num_genes)){
  stop("--num-genes must be an integer.")
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
library_ids_sce_order <- stringr::str_extract(sce_files, pattern = library_search)

sce_file_df <- data.frame(sce_files = sce_files, 
                          library_id= library_ids_sce_order) %>%
  dplyr::left_join(library_metadata_df, by = c("library_id" = "library_biomaterial_id")) %>%
  dplyr::select(library_id, opt$grouping_var, sce_files)

# group dataframe by the grouping variable 
grouped_sce_file_df <- split(sce_file_df, sce_file_df[,opt$grouping_var])

# create a list of SCE lists that is named by the grouping variable with 
# each individual inner SCE list named by the library IDs
create_grouped_sce_list <- function(sce_info_dataframe){
  
  library_sce_list = list()
  for (library_idx in 1:length(sce_info_dataframe$library_id)){
    
    # read sce list for each library
    sce <- readr::read_rds(sce_info_dataframe$sce_files[library_idx])
    library_name <- sce_info_dataframe$library_id[library_idx]
    # create a list for each group named by the library IDs
    library_sce_list[[library_name]] <- sce
    
  }
  
  return(library_sce_list)
  
}  

grouped_sce_list <- grouped_sce_file_df %>%
  purrr::map(create_grouped_sce_list)
  

# create a list of merged SCE objects by group
merged_sce_list <- grouped_sce_list %>%
  purrr::map(~ combine_sce_objects(.x, preserve_rowdata_columns = c("Gene", "gene_names")))

# HVG and dim reduction --------------------------------------------------------

# apply HVG calculation to list of merged SCEs
# object will only be subset to HVG if subset_hvg is true
merged_sce_list <- merged_sce_list %>%
  purrr::map(~ add_var_genes(.x, 
                             num_genes = opt$num_genes,
                             subset_hvg = opt$subset_hvg))

# add PCA and UMAP 
# if --use_all_genes is used, use all genes, otherwise only HVG are used
if(opt$use_all_genes){
  merged_sce_list <- merged_sce_list %>%
    purrr::map( ~ perform_dim_reduction(.x, 
                                        var_genes = rownames(.x),
                                        pca_type = "multi"))
} else {
  merged_sce_list <- merged_sce_list %>%
    purrr::map( ~ perform_dim_reduction(.x, 
                                        var_genes = metadata(.x)$variable_genes,
                                        pca_type = "multi")) 
}


# Write RDS --------------------------------------------------------------------

# create paths to merged SCE files
# named with the name of the sce list which corresponds to the grouping variable, not library ID
merged_sce_files <- file.path(opt$merged_sce_dir,
                              paste0(names(merged_sce_list), 
                                     "_merged_sce.rds"))

# export files 
purrr::walk2(merged_sce_list, merged_sce_files, readr::write_rds)
