# Script used to perform integration on SCE objects using R-based methods
#
# Merged SCE files, which have been grouped by a specified grouping variable that 
# must be present as a column in the `library_file`, are read in. Integration is 
# performed with either `fastMNN` or `harmony`, depending on user specification, 
# and integrated SCE objects are stored as RDS files in the `--integrated_sce_dir`.

# Option descriptions: 
#
# --library_file: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column
# --grouping_var: Column name present in the library metadata file that identifies 
#   combined SCE objects
# --merged_sce_dir: Path to folder where all merged SCE objects have been stored 
# --integrated_sce_dir: Path to folder where all integrated SCE objects will be stored. 
#   Inside this folder are nested folders for results from each integration method
# --batch_column: Column name present in SCE object colData slot that indicates batch
#   groupings
# --method: The integration method to use, either `fastMNN` or `harmony` (case-insensitive).
# --seed: Random seed to use during integration
# --harmony_covariate_cols: Optional comma-separated list of additional covariate columns to consider
#   during integration with `harmony`. This is argument is ignored if the provided
#   method is `fastMNN`
# --harmony_from_expression: Flag to specify that `harmony` integration should be
#   performed from the expression matrix directly rather than from pre-computed PCs.
#   This is argument is ignored if the provided method is `fastMNN`
#   This is argument is ignored if the provided method is `fastMNN`
# --fastmnn_no_cosine: Flag to specify that `fastMNN` integration should not perform
#   cosine normalization. This is argument is ignored if the provided method is `harmony`
# --fastmnn_gene_list: Optional comma-separated list of genes for `fastMNN` to consider
#    during integration, if all genes are not desired.
# --integration_options: CURRENTLY NOT USED. Additional options to pass into the specified
#   method 
# 
# Usage examples:
# Rscript 03-integrate-sce.R --method=fastMNN   # Default fastMNN usage
# Rscript 03-integrate-sce.R --method=fastMNN --fastmnn_no_cosine   # fastMNN without cosine normalization
#
# Rscript 03-integrate-sce.R --method=harmony   # Default harmony usage
# Rscript 03-integrate-sce.R --method=harmony --harmony_from_expression   # harmony starting from expression matrix




# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
library(magrittr)
library(optparse)

# source helper functions and integration functions
source(file.path(project_root, "scripts", "utils", "integration-helpers.R"))
source(file.path(project_root, "scripts", "utils", "integrate-fastMNN.R"))
source(file.path(project_root, "scripts", "utils", "integrate-harmony.R"))

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
    help = "Column name present in the library metadata file that identifies combined SCE objects"
  ),
  make_option(
    opt_str = c("--merged_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "merged-sce-objects"),
    help = "path to folder where all merged SCE objects are stored as RDS files"
  ),
  make_option(
    opt_str = c("--integrated_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "human_cell_atlas", "integrated-sce-objects"),
    help = "path to folder where all integrated SCE objects will be saved as RDS files"
  ),
  make_option(
    opt_str = c("--method"),
    type = "character",
    default = NULL,
    help = "Integration method to use, either `fastMNN` or `harmony` (case-insensitive)."
  ),
  make_option(
    opt_str = c("-b", "--batch_column"), 
    type = "character",
    default = "batch", 
    help = "Name of the column in the SCE object that indicates batch groupings."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2022,
    help = "random seed to set during integration"
  ),
  make_option(
    opt_str = c("--harmony_covariate_cols"), 
    type = "character",
    default = NULL, 
    help = "Optional comma-separated list of columns (e.g. patient, sex) to consider as covariates
            during integration with `harmony`."
  ),
  make_option(
    opt_str = c("--harmony_from_expression"),
    action = "store_false",
    help = "Indicate whether to use the gene expression matrix, rather than PCs, during `harmony` integration.
    To use expression instead of PCs (default), use `--harmony_from_expression`"
  ),
  make_option(
    opt_str = c("--fastmnn_no_cosine"),
    action = "store_false",
    help = "Indicate whether to turn off cosine normalization during `fastMNN` integration. 
    To turn off cosine normalization, use `--fastmnn_no_norm`"
  ),
  make_option(
    opt_str = c("--fastmnn_gene_list"),
    type = "character",
    default = NULL,
    help = "Optional comma-separated list of genes for `fastMNN` to consider during integration, 
    if all genes are not desired."
  ),
  make_option(
    opt_str = c("--integration_options"),
    type = "character",
    default = NULL,
    help = "CURRENTLY NOT USED. Additional options to pass into the given integration method"
  )
)

# Setup ------------------------------------------------------------------------
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# checks that provided metadata files exist
if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# read in library metadata to find `grouping_var` names, which corresponds to
#  the name of saved merged objects as: `{grouping_var}_merged_sce.rds`
library_metadata_df <- readr::read_tsv(opt$library_file)
grouping_var <- as.symbol(opt$grouping_var)
group_names <- library_metadata_df %>%
  # arrange to ensure proper mapping when naming merged list later
  dplyr::arrange({{grouping_var}}) %>%
  dplyr::pull({{grouping_var}}) %>%
  unique()



# Check and assign provided method
if (is.null(opt$method)) {
  stop("You must specify either `fastMNN` or `harmony` (case-insensitive) to --method.")
} else {
  integration_method <- tolower(opt$method)
  if (!(integration_method %in% c("fastmnn", "harmony"))) {
    stop("You must specify either `fastMNN` or `harmony` (case-insensitive) to --method.")
  }
}


# setup output directory, which is opt$integrated_sce_dir and a subdirectory for 
#  the given method
if(!dir.exists(opt$integrated_sce_dir)){
  dir.create(opt$integrated_sce_dir, recursive = TRUE)
}
method_integrated_sce_dir <- file.path(opt$integrated_sce_dir, integration_method)
if(!dir.exists(method_integrated_sce_dir)){
  dir.create(method_integrated_sce_dir, recursive = TRUE)
}

# Identify SCE files -----------------------------------------------------------

# find SCE files that match library ID 
# the only files in this directory should be RDS
group_search <- paste(group_names, collapse = "|")
merged_sce_files <- list.files(opt$merged_sce_dir,
                        pattern = group_search, 
                        recursive = TRUE,
                        full.names = TRUE)

# if the number of sce files is different than the project names, find the missing files
if(length(merged_sce_files) < length(group_names)){
  
  groups_found <- stringr::str_extract(merged_sce_files, group_search)
  missing_groups <- setdiff(group_names, groups_found)
  
  stop(
    glue::glue(
      "\nMissing merged SCE object(s) for {missing_projects}.
      Make sure that you have run `02-prepare-merged-sce.R`."
    )
  )
}


# Read in all merged SCE files -----------------------------
merged_sce_objs <- merged_sce_files %>%
  purrr::map(readr::read_rds)


# Define helper function for splitting up comma-separated arguments ------

split_comma_args <- function(arg) {
  arg %>%
    # Remove spaces
    stringr::str_replace_all(" ", "") %>%
    # Split on ,
    stringr::str_split(",") %>%
    # magrittr for the win!
    extract2(1)
}

# Perform integration with fastMNN, if specified -------------------------

if (integration_method == "fastmnn") {

  # Set up `cosine_norm` argument based on user options
  if (is.null(opt$fastmnn_no_norm)) {
    fastmnn_cosine_norm <- TRUE
  } else {
    fastmnn_cosine_norm <- FALSE
  }
  # Set up `gene_list` argument based on user options
  if (is.null(opt$fastmnn_gene_list)) {
    fastmnn_gene_list <- NULL
  } else {
    fastmnn_gene_list <- split_comma_args(opt$fastmnn_gene_list)
  }
    
    
  # Run fastMNN, with or without additional options passed in
 # if (is.null(opt$integration_options)) {
    integrated_sce_objs <- merged_sce_objs %>%
      purrr::map(
        ~integrate_fastMNN(
          .x,
          batch_column = opt$batch_column,
          cosine_norm  = fastmnn_cosine_norm,
          gene_list    = fastmnn_gene_list,
          seed         = opt$seed
        ))
 # } else {
 #   # DOES NOT WORK, HOW TO PASS IN `...` FROM OPTPARSE ARGS?
 #   integrated_sce_objs <- merged_sce_objs %>%
 #     purrr::map(
 #       ~integrate_fastMNN(
 #         .x,
 #         batch_column = opt$batch_column,
 #         cosine_norm  = fastmnn_cosine_norm,
 #         gene_list    = fastmnn_gene_list,
 #         seed         = opt$seed, 
 #         opt$integration_options
 #       ))
 # }
}

# Perform integration with harmony, if specified -------------------------

if (integration_method == "harmony") {

  
  # Set up `from_pca` argument based on user options
  if (is.null(opt$harmony_from_expression)) {
    harmony_from_pca <- TRUE
  } else {
    harmony_from_pca <- FALSE
  }
  
  # Set up `covariate_cols` argument based on user options
  if (is.null(opt$harmony_covariate_cols)) {
    harmony_covariate_cols <- c()
  } else {
    harmony_covariate_cols <- split_comma_args(opt$harmony_covariate_cols)
  }

  # Run harmony, with or without additional options passed in
  #if (is.null(opt$integration_options)) {
    integrated_sce_objs <- merged_sce_objs %>%
      purrr::map(
        ~integrate_harmony(
          .x,
          batch_column = opt$batch_column,
          covariate_cols = harmony_covariate_cols,
          from_pca = harmony_from_pca,
          seed = opt$seed
      ))
  #} else {
  #  # DOES NOT WORK, HOW TO PASS IN `...` FROM OPTPARSE ARGS?
  #  integrated_sce_objs <- merged_sce_objs %>%
  #    purrr::map(
  #      ~integrate_harmony(
  #        .x,
  #        batch_column = opt$batch_column,
  #        covariate_cols = harmony_covariate_cols,
  #        from_pca = harmony_from_pca,
  #        seed = opt$seed,
  #        opt$integration_options
  #      ))
  #}
  
}


# Write integrated SCE objects to RDS --------------------------------------------------------------------

# create paths to integrated SCE files
integrated_sce_files <- file.path(method_integrated_sce_dir,
                              paste0(group_names, 
                                     "_integrated_sce.rds"))

# export files 
purrr::walk2(integrated_sce_objs, integrated_sce_files, readr::write_rds)
