# Script used to prepare SCE objects to be integrated
#
# SCE files are grouped by a specified grouping variable that must be present as
# a column in the `library_file`. SCE objects for each group are merged, highly
# variable genes for the merged object are calculated, and then merged objects are
# stored as RDS files in the `--merged_sce_dir`. This script also adds a column
# with cell type information to each SCE object prior to merging, if cell type 
# information is available.

# Option descriptions:
#
# --library_file: The path to the file listing all libraries that should be included
#   in the conversion from loom to SCE. This file must contain the
#   `library_biomaterial_id` column
# --grouping_var: Column name present in the library metadata file to use for
#   grouping SCE objects and merging.
# --add_celltype: Boolean indicating whether or not celltype data should be added to the 
#   individual SCE object prior to merging. 
# --celltype_info: Path to file containing cell type information for each SCE object. 
#   Must contain columns named `library_biomaterial_id`, `celltype`, and `barcode`.
#   Only required if --add_celltype is used.
# --batch_column: The name of the column in colData that indicates the batches for each cell,
#   typically this corresponds to the library id. Default is "batch".
# --random_merge: Used to indicate whether or not to merge SCE objects in a random order. 
#   Default is FALSE.
# --subset_hvg: Indicates whether or not to subset the merged SCE object by highly variable genes.
#   If --subset_hvg is used, the merged SCE object will only contain genes
#   identified as highly variable genes.
# --pca_use_all_genes: Indicates whether or not to use the all genes as input to performing
#   principal component analysis. Otherwise only highly variable genes are used
#   as input.
# --num_hvg: Number of highly variable genes to select.
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
    opt_str = c("--groups_to_integrate"),
    type = "character",
    default = "All",
    help = "Groups present in `grouping_var` column of metadata file to create merged SCE objects for.
      Default is 'All' which produces a merged object for each group in the metadata file. 
      Specify groups by using a vector, e.g. c('group1','group2')"
  ),
  make_option(
    opt_str = c("--add_celltype"),
    action = "store_true",
    default = FALSE,
    help = "Boolean indicating whether or not celltype data should be added to the 
     individual SCE object prior to merging."
  ),
  make_option(
    opt_str = c("--celltype_info"),
    type = "character",
    default = file.path(project_root, "sample-info", "hca-celltype-info.tsv"),
    help = "Path to file containing cell type information for each SCE object. 
      Must contain columns named `library_biomaterial_id`, `celltype`, and `barcode`."
  ),
  make_option(
    opt_str = c("--batch_column"),
    type = "character",
    default = "batch",
    help = "The name of the column in colData that indicates the batches for each cell,
      typically this corresponds to the library id. Default is 'batch'."
  ),
  make_option(
    opt_str = c("--random_merge"),
    default = FALSE,
    action = "store_true",
    help = "Used to indicate whether or not to merge SCE objects in a random order. Default is FALSE."
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
    opt_str = c("--pca_use_all_genes"),
    default = FALSE,
    action = "store_true",
    help = "Indicates whether or not to use the all genes as input to performing
      principal component analysis. Otherwise only highly variable genes are used
      as input."
  ),
  make_option(
    opt_str = c("-n", "--num_hvg"),
    dest = "num_genes",
    type = "integer",
    default = 5000,
    help = "Number of highly variable genes to select."
  ),
  make_option(
    opt_str = c("--merged_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "merged_sce"),
    help = "path to folder where all merged SCE objects files will be saved as RDS files"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = NULL,
    help = "random seed to set prior to merging"
  )
)

# Setup ------------------------------------------------------------------------

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)

# check that num genes provided is an integer
if(!is.integer(opt$num_genes)){
  stop("--num_hvg must be an integer.")
}

# checks that provided metadata files exist
if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# read in library metadata and grab unfiltered sce file paths
library_metadata_df <- readr::read_tsv(opt$library_file)

# check that cell type file exists if using add_celltype option 
if(opt$add_celltype){
  if(!file.exists(opt$celltype_info)){
    stop("--celltype_info file provided does not exist.")
  }
  celltype_info_df <- readr::read_tsv(opt$celltype_info)
}

# check that grouping variable is present
if(!opt$grouping_var %in% colnames(library_metadata_df)){
  stop("Must provide a grouping_var that is a column in the library metadata file.")
}

# define groups to integrate
groups <- library_metadata_df %>%
  dplyr::pull(opt$grouping_var) %>%
  unique()

# check that groups to integrate isn't set to All 
if(length(opt$groups_to_integrate) == 1 && (opt$groups_to_integrate == "All")){
  groups_to_integrate <- groups
} else {
  # if All is not used then unlist groups, using space, needed to parse a list from snakemake
  groups_to_integrate <- unlist(stringr::str_split(opt$groups_to_integrate, " "))
  
  # check that specified groups are present in grouping_var column 
  if(!any(groups_to_integrate %in% groups)){
    stop("Provided `--groups_to_integrate` must also be present in the `--grouping_var` column of 
         the library metadata file.")
  }
}

# subset metadata file to only contain groups to integrate 
library_metadata_df <- library_metadata_df %>%
  dplyr::filter(.data[[opt$grouping_var]] %in% groups_to_integrate)

# grab library ids 
library_ids <- library_metadata_df %>%
  dplyr::pull(library_biomaterial_id)

# setup output directory
create_dir(opt$merged_sce_dir)

# Identify SCE files -----------------------------------------------------------

# grab all unique directories corresponding to projects considered for integration
input_sce_dirs <- library_metadata_df %>%
  dplyr::pull(integration_input_dir) %>%
  unique()

# find SCE files that match library ID, and throw an error if any are missing.
sce_files <- find_input_sce_files(library_ids, input_sce_dirs)

# Merge by group ---------------------------------------------------------------

# get the library IDs from the SCE file names so that we can name the SCEs in the correct order
library_ids_sce_order <- stringr::str_extract(sce_files, 
                                              pattern = paste(library_ids, collapse = "|"))

sce_file_df <- data.frame(sce_files = sce_files,
                          library_id= library_ids_sce_order) %>%
  dplyr::left_join(library_metadata_df, by = c("library_id" = "library_biomaterial_id")) %>%
  dplyr::select(library_id, opt$grouping_var, sce_files)

# group dataframe by the grouping variable
grouped_sce_file_df <- split(sce_file_df, sce_file_df[,opt$grouping_var])

# create a list of SCE lists that is named by the grouping variable with
# each individual inner SCE list named by the library IDs
create_grouped_sce_list <- function(sce_info_dataframe,
                                    celltype_info_df = NULL,
                                    add_celltype){

  library_sce_list = list()
  for (library_idx in 1:length(sce_info_dataframe$library_id)){

    # read sce list for each library
    sce <- readr::read_rds(sce_info_dataframe$sce_files[library_idx])
    library_name <- sce_info_dataframe$library_id[library_idx]
    
    if(add_celltype){
      # if celltype info is provided add to sce object 
      if(!is.null(celltype_info_df)){
        # check that library has cell type information
        if(library_name %in% unique(celltype_info_df$library_biomaterial_id)){
          
          # filter celltype info to only have info for specified library 
          filtered_celltype_info <- celltype_info_df %>%
            dplyr::filter(library_biomaterial_id == library_name) %>%
            dplyr::select(barcode, celltype)
          
          # add celltype info 
          sce <- add_celltype_info(sce_object = sce, 
                                   celltype_info_df = filtered_celltype_info) 
          
          # add flag indicating that cell type information is available 
          metadata(sce)$celltype_info_available <- TRUE 
        }
        
        # only add celltype column/metadata if add celltype is yes, but no celltype data is available 
      } else {
        # note that no cell type information is available
        colData(sce)$celltype <- NA
        metadata(sce)$celltype_info_available <- FALSE
      } 
    }
    
    # create a list for each group named by the library IDs
    library_sce_list[[library_name]] <- sce

  }

  return(library_sce_list)

}

# create grouped sce list with/without celltype addition
if(opt$add_celltype){
  grouped_sce_list <- grouped_sce_file_df %>%
    purrr::map(create_grouped_sce_list, celltype_info_df, add_celltype = opt$add_celltype) 
} else {
  grouped_sce_list <- grouped_sce_file_df %>%
    purrr::map(create_grouped_sce_list, add_celltype = opt$add_celltype)
}

# create a list of merged SCE objects by group
#  In this default usage, a batch column named `batch` will get created
merged_sce_list <- grouped_sce_list %>%
  purrr::map(combine_sce_objects, 
             preserve_rowdata_columns = c("Gene", "gene_names", "ensembl_ids"),
             random_merge = opt$random_merge)

# HVG and dim reduction --------------------------------------------------------

# apply HVG calculation to list of merged SCEs
# object will only be subset to HVG if subset_hvg is true
merged_sce_list <- merged_sce_list %>%
  purrr::map(~ set_var_genes(.x,
                             num_genes = opt$num_genes,
                             subset_hvg = opt$subset_hvg,
                             batch_column = opt$batch_column))

# add PCA and UMAP
# if --pca_use_all_genes is used, use all genes, otherwise only HVG are used
if(opt$pca_use_all_genes){
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
