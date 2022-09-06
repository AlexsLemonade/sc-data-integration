# Script used to perform integration on a given merged SCE object using R-based methods
#
# A merged SCE file in RDS format is read in. Integration is performed with either
# `fastMNN`, `harmony`, or `Seurat` (using either `cca` or `rpca` reduction)
# depending on user specification. The integrated SCE object is saved as an RDS
# file in the provided output file.
#
# Option descriptions:
#
# --input_sce_file: Path to RDS file with merged SCE object to integrate
# --output_sce_file: Path to RDS file where the integrated SCE object will be saved.
#   Inside this folder are nested folders for results from each integration method
# --method: The integration method to use, either `fastMNN` or `harmony` (case-insensitive).
# --batch_column: Column name present in SCE object colData slot that indicates batch
#   groupings
# --seed: Random seed to use during integration
# --harmony_covariate_cols: Optional comma-separated list of additional covariate columns to consider
#   during integration with `harmony`. This argument is ignored if the provided
#   method is not `harmony`
# --harmony_from_expression: Flag to specify that `harmony` integration should be
#   performed from the expression matrix directly rather than from pre-computed PCs.
#   This argument is ignored if the provided method is not `harmony`
# --fastmnn_no_cosine: Flag to specify that `fastMNN` integration should not perform
#   cosine normalization. This argument is ignored if the provided method is not
#  `fastMNN`
# --fastmnn_use_all_genes: Flag to specify that `fastMNN` integration should all genes in
#   during integration instead of, by default, using only the previously-identified HVGs,
#   stored in the `variable_genes` column of metadata slot. This argument is ignored if
#   the provided method is not `fastMNN`.
# --seurat_reduction_method: Seurat reduction method to use, either `cca` or `rcpa`.
#   This argument is ignored if the provided method is not `seurat`. There is no 
#   default, so this argument is required if the provided method is `seurat`.
# --seurat_num_genes: Number of variables genes Seurat should identify.
#   This argument is ignored if the provided method is not `seurat`. Default: 2000.
# --seurat_integration_dims: Number of dimensions Seurat should use during integration.
#   This argument is ignored if the provided method is not `seurat`. Default:30.
# --seurat_umap_dims: Number of dimensions Seurat should use during UMAP.
#   This argument is ignored if the provided method is not `seurat`. Default:30.
# --corrected_only: Flag to specify that only corrected gene expression values should
#   be returned in the integrated SCE object. Default usage of this script will
#   return all data.
#
#
# Usage examples:
#
# Rscript 03-integrate-sce.R \
#  --input_sce_file ../results/human_cell_atlas/merged-sce/1M_Immune_Cells_merged_sce.rds \
#  --output_sce_file ../results/human_cell_atlas/integrated_sce/1M_Immune_Cells_integrated_harmony_sce.rds \
#  --method=harmony
#
# Rscript 03-integrate-sce.R \
#  --input_sce_file ../results/human_cell_atlas/merged-sce/1M_Immune_Cells_merged_sce.rds \
#  --output_sce_file ../results/human_cell_atlas/integrated-sce/fastMNN/1M_Immune_Cells_integrated_fastmnn_sce.rds \
#  --method=fastmnn
#
# Rscript 03-integrate-sce.R \
#  --input_sce_file ../results/human_cell_atlas/merged-sce/1M_Immune_Cells_merged_sce.rds \
#  --output_sce_file ../results/human_cell_atlas/integrated-sce/fastMNN/1M_Immune_Cells_integrated_seurat-rpca_sce.rds \
#  --method=seurat --seurat_reduction_method=rpca
#
# load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
library(magrittr)
library(optparse)

# source integration functions
utils_path <- file.path(project_root, "scripts", "utils")
suppressPackageStartupMessages({
  source(file.path(utils_path, "integrate-fastMNN.R"))
  source(file.path(utils_path, "integrate-harmony.R"))
  source(file.path(utils_path, "integrate-seurat.R"))
  source(file.path(utils_path, "integration-helpers.R"))
})

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--input_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file that contains the merged SCE object to integrate"
  ),
  make_option(
    opt_str = c("--output_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file where the integrated SCE object will be saved"
  ),
  make_option(
    opt_str = c("--method"),
    type = "character",
    default = NULL,
    help = "Integration method to use, either `fastMNN`, `harmony`, or `Seurat` (case-insensitive)."
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
    opt_str = c("--seurat_reduction_method"),
    type = "character",
    default = NULL,
    help = "Reduction method to use during Seurat integration."
  ),
  make_option(
    opt_str = c("--seurat_num_genes"),
    type = "numeric",
    default = 2000,
    help = "Number of variable genes for Seurat to identify and use."
  ),
  make_option(
    opt_str = c("--seurat_integration_dims"),
    type = "numeric",
    default = 30,
    help = "Number of dimensions for Seurat to use during integration."
  ),
  make_option(
    opt_str = c("--seurat_umap_dims"),
    type = "numeric",
    default = 30,
    help = "Number of dimensions for Seurat to use during UMAP calculation."
  ),
  make_option(
    opt_str = c("--corrected_only"),
    action = "store_true",
    default = FALSE,
    help = "Indicate whether to only return the corrected gene expression data, and
    not uncorrected expression, in the integrated SCE object. To return only
    corrected expression, use `--corrected_only`."
  )
)


# --corrected_only: Flag to specify that only corrected gene expression values should
#   be returned in the integrated SCE object. Default usage of this script will
#   return all data.

# Setup ------------------------------------------------------------------------
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))



# Check and assign provided method based on available methods
available_methods <- c("fastmnn", "harmony", "seurat")

# helper function for method check fails
stop_method <- function() {
  stop(
    paste("You must specify one of the following (case-insensitive) to --method:",
          paste(available_methods, collapse = ", ")
    )
  )
}

if (is.null(opt$method)) {
  stop_method()
} else {
  integration_method <- tolower(opt$method)
  if (!(integration_method %in% available_methods)) {
    stop_method()
  }
}


# Check that provided input file exists and is an RDS file
if(is.null(opt$input_sce_file)) {
  stop("You must provide the path to the RDS file with merged SCEs to --input_sce_file")
} else {
  if(!file.exists(opt$input_sce_file)) {
    stop("Provided --input_sce_file file does not exist. Make sure that you have
          run `02-prepare-merged-sce.R` or the provided path is correct.")
  }
}

# Check that both input and output files have RDS extensions
if(!(grepl("\\.rds$", opt$input_sce_file, ignore.case = TRUE)) ||
   !(grepl("\\.rds$", opt$output_sce_file, ignore.case = TRUE))) {
  stop("The provided --input_sce_file and --output_sce_file files must be RDS files.")
}


# Check that directory for output file exists and the specified file is an RDS file
integrated_sce_dir <- dirname(opt$output_sce_file)
if(!dir.exists(integrated_sce_dir)){
  dir.create(integrated_sce_dir, recursive = TRUE)
}



# Read in SCE file -----------------------------
merged_sce_obj <- readr::read_rds(opt$input_sce_file)


# Perform integration with fastMNN, if specified -------------------------
if (integration_method == "fastmnn") {

  # Prepare `gene_list` argument
  if (opt$fastmnn_use_all_genes) {
    fastmnn_gene_list <- NULL
  } else {
    fastmnn_gene_list <- metadata(merged_sce_obj)$variable_genes
  }

  # Perform integration
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

  # Perform integration
  integrated_sce_obj <- integrate_harmony(
    merged_sce_obj,
    batch_column = opt$batch_column,
    covariate_cols = opt$harmony_covariate_cols,
    from_pca = opt$harmony_from_expression,
    seed = opt$seed
  )
}

# Perform integration with seurat, if specified -------------------------
if (integration_method == "seurat") {

  # Convert the merged sce object into a list of seurat objects
  seurat_list <- prepare_seurat_list(merged_sce_obj)
  
  # Perform integration
  integrated_seurat_obj <- integrate_seurat(
    seurat_list,
    opt$seurat_reduction_method, # integrate_seurat() will check this argument
    batch_column = opt$batch_column,
    num_genes = opt$seurat_num_genes,
    integration_dims = 1:opt$seurat_integration_dims,
    umap_dims = 1:opt$seurat_umap_dims
  )
  
  # Converted the integrated seurat object into an SCE object
  integrated_sce_obj <- as.SingleCellExperiment(integrated_seurat_obj)
  
  # Rename assays appropriately:
  
  # The converted `logcounts` is the CORRECTED expression
  assay(integrated_sce_obj, paste0(opt$seurat_reduction_method, "_corrected")) <- logcounts(integrated_sce_obj)
  
  # The logcounts in SCT assay is logcounts
  logcounts(integrated_sce_obj) <- logcounts( altExp(integrated_sce_obj, "SCT") )
  
  # The counts in RNA assay is counts  
  counts(integrated_sce_obj) <- counts( altExp(integrated_sce_obj, "RNA") )
  
  # Remove altExps
  integrated_sce_obj <- removeAltExps(integrated_sce_obj)

  # Convert reducedDims names back to <lowercase>_<UPPERCASE> because 
  #   `as.SingleCellExperiment` makes them all uppercase
  reducedDimNames(integrated_sce_obj) <- stringr::str_replace_all(
    reducedDimNames(integrated_sce_obj),
    toupper(opt$seurat_reduction_method),
    opt$seurat_reduction_method
  )

}


# Remove uncorrected expression values, if specified ----------
if (opt$corrected_only) {
  integrated_sce_obj <- remove_uncorrected_expression(integrated_sce_obj)
}


# Write integrated SCE object to RDS -------------------
readr::write_rds(integrated_sce_obj, opt$output_sce_file)
