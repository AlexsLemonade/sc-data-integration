# Script used to perform integration on a given merged SCE object using R-based methods
#
# A merged SCE file, which has been grouped by a specified grouping variable that 
# must be present as a column in the `library_file`, is read in. Integration is 
# performed with either `fastMNN` or `harmony`, depending on user specification, 
# and the integrated SCE object is saved as RDS files in a nested folder within
# `--integrated_sce_dir` named for the given integration method.

# Option descriptions: 
#
# --merged_sce_file: Path to RDS file with merged SCE object to integrate
# --integrated_sce_file: Path to RDS file where the integrated SCE object will be saved. 
#   Inside this folder are nested folders for results from each integration method
# --method: The integration method to use, either `fastMNN` or `harmony` (case-insensitive).
# --batch_column: Column name present in SCE object colData slot that indicates batch
#   groupings
# --seed: Random seed to use during integration
# --harmony_covariate_cols: Optional comma-separated list of additional covariate columns to consider
#   during integration with `harmony`. This argument is ignored if the provided
#   method is `fastMNN`
# --harmony_from_expression: Flag to specify that `harmony` integration should be
#   performed from the expression matrix directly rather than from pre-computed PCs.
#   This argument is ignored if the provided method is `fastMNN`
# --fastmnn_no_cosine: Flag to specify that `fastMNN` integration should not perform
#   cosine normalization. This argument is ignored if the provided method is `harmony`
# --fastmnn_use_all_genes: Flag to specify that `fastMNN` integration should all genes in
#   during integration instead of, by default, using only the previously-identified HVGs, 
#   stored in the `variable_genes` column of metadata slot. This argument is ignored if 
#   the provided method is `harmony`
# 
# 
# Usage examples:
#
# Rscript 03-integrate-sce.R \
#  --merged_sce_file ../results/human_cell_atlas/merged-sce-objects/1M_Immune_Cells_merged_sce.rds \
#  --integrated_sce_file ../results/human_cell_atlas/integrated-sce-objects/harmony/1M_Immune_Cells_merged_sce.rds \
#  --method=harmony 
#
# Rscript 03-integrate-sce.R \
#  --merged_sce_file ../results/human_cell_atlas/merged-sce-objects/1M_Immune_Cells_merged_sce.rds \
#  --integrated_sce_file ../results/human_cell_atlas/integrated-sce-objects/fastMNN/1M_Immune_Cells_merged_sce.rds \
#  --method=fastmnn
#

# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
library(magrittr)
library(optparse)

# source integration functions
source(file.path(project_root, "scripts", "utils", "integrate-fastMNN.R"))
source(file.path(project_root, "scripts", "utils", "integrate-harmony.R"))

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--merged_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file that contains the merged SCE object to integrate"
  ),
  make_option(
    opt_str = c("--integrated_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file where the integrated SCE object will be saved"
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
    default = NULL,
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
    default = TRUE,
    help = "Indicate whether to use the gene expression matrix, rather than PCs, during `harmony` integration.
    To use expression instead of PCs (default), use `--harmony_from_expression`"
  ),
  make_option(
    opt_str = c("--fastmnn_no_cosine"),
    action = "store_false",
    default = TRUE,
    help = "Indicate whether to turn off cosine normalization during `fastMNN` integration. 
    To turn off cosine normalization, use `--fastmnn_no_cosine`"
  ),
  make_option(
    opt_str = c("--fastmnn_use_all_genes"),
    action = "store_true", 
    default = FALSE,
    help = "Indicate whether to use all genes instead of only HVGs during `fastMNN` integration. 
    To use all genes, use `--fastmnn_use_all_genes`"
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

# Check and assign provided method
if (is.null(opt$method)) {
  stop("You must specify either `fastMNN` or `harmony` (case-insensitive) to --method.")
} else {
  integration_method <- tolower(opt$method)
  if (!(integration_method %in% c("fastmnn", "harmony"))) {
    stop("You must specify either `fastMNN` or `harmony` (case-insensitive) to --method.")
  }
}

# Check that provided input file exists
if(is.null(opt$merged_sce_file)) {
  stop("You must provide the path to the RDS file with merged SCEs to --merged_sce_file")
} else {
  if(!file.exists(opt$merged_sce_file)) {
    stop("Provided --merged_sce_file file does not exist. Make sure that you have 
          run `02-prepare-merged-sce.R` or the provided path is correct.")
  }
}

# Check that directory for output file exists
integrated_sce_dir <- dirname(opt$integrated_sce_file)
if(!dir.exists(integrated_sce_dir)){
  dir.create(integrated_sce_dir, recursive = TRUE)
}

# Read in SCE file -----------------------------
merged_sce_obj <- readr::read_rds(opt$merged_sce_file)


# Perform integration with fastMNN, if specified -------------------------

if (integration_method == "fastmnn") {
    
  # Prepare `gene_list` argument
  if (opt$fastmnn_use_all_genes) {
    fastmnn_gene_list <- NULL
  } else {
    fastmnn_gene_list <- metadata(merged_sce_obj)$variable_genes
  }
  
  # Run fastMNN, without additional options passed in
  integrated_sce_obj <- integrate_fastMNN(
    merged_sce_obj,
    batch_column = opt$batch_column,
    cosine_norm  = opt$fastmnn_no_cosine,
    gene_list    = fastmnn_gene_list,
    seed         = opt$seed
  )

}

# Perform integration with harmony, if specified -------------------------

if (integration_method == "harmony") {

  
  # Set up `covariate_cols` argument based on user options
  if (is.null(opt$harmony_covariate_cols)) {
    harmony_covariate_cols <- c()
  } else {
    harmony_covariate_cols <- unlist(stringr::str_split(opt$harmony_covariate_cols, "\\s*,\\s*"))
  }

  # Run harmony, without additional options passed in
  #if (is.null(opt$integration_options)) {
    integrated_sce_obj <- integrate_harmony(
      merged_sce_obj,
      batch_column = opt$batch_column,
      covariate_cols = opt$harmony_covariate_cols,
      from_pca = opt$harmony_from_expression,
      seed = opt$seed
    )
  #} else {
  #  # DOES NOT WORK
  #  integrated_sce_obj <- integrate_harmony(
  #    merged_sce_obj,
  #    batch_column = opt$batch_column,
  #    covariate_cols = harmony_covariate_cols,
  #    from_pca = harmony_from_pca,
  #    seed = opt$seed,
  #    opt$integration_options
  #  )
  #}
  
}


# Write integrated SCE object to RDS -------------------

readr::write_rds(integrated_sce_obj, opt$integrated_sce_file)
