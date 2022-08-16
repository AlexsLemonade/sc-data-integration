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
    opt_str = c("--input_h5_file"),
    type = "character",
    default = "results/human_cell_atlas/integrated_anndata/1M_Immune_Cells_integrated_scVI.h5",
    help = "Path to H5 file that contains the integrated AnnData object"
  ),
  make_option(
    opt_str = c("--output_sce_file"),
    type = "character",
    default = "results/human_cell_atlas/integrated_sce/1M_Immune_Cells_integrated_scVI_sce.rds",
    help = "Path to RDS file where the integrated SCE object will be saved"
  ),
  make_option(
    opt_str = c("--method"),
    type = "character",
    default = "scVI",
    help = "Method that was used for integration, either `scanorama` or `scVI` (case-insensitive)."
  ),
  make_option(
    opt_str = c("-b", "--batch_column"),
    type = "character",
    default = "batch",
    help = "Name of the column in the AnnData object that indicates batch groupings."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = NULL,
    help = "random seed to set"
  ),
  make_option(
    opt_str = c("--corrected_only"),
    action = "store_true",
    default = FALSE,
    help = "Indicate whether the integrated AnnData object contains only 
    corrected gene expression data, and is missing uncorrected expression (no `X` matrix). 
    To indicate the object only contains corrected data, use `--corrected_only`."
  )
)

# Setup ------------------------------------------------------------------------
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Check provided method based on available methods
available_methods <- c("scanorama", "scvi")

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
if(is.null(opt$input_h5_file)) {
  stop("You must provide the path to the H5 file with merged SCEs to --input_h5_file")
} else {
  if(!file.exists(opt$input_h5_file)) {
    stop("Provided --input_h5_file file does not exist.")
  }
}

# Check that both input and output files have correct extensions
if(!(stringr::str_ends(opt$input_h5_file, ".hdf5|.h5"))){
  stop("`--input_h5_file` must end in either '.hdf5' or '.h5'")
}
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("`--output_sce_file` must end in '.rds'")
}

# Check that directory for output file exists
integrated_sce_dir <- dirname(opt$output_sce_file)
if(!dir.exists(integrated_sce_dir)){
  dir.create(integrated_sce_dir, recursive = TRUE)
}

# Anndata to SCE post-processing -----------------------------------------------

# read in AnnData object as H5
integrated_sce <- zellkonverter::readH5AD(opt$input_h5_file, 
                                          # if corrected only is set, skips reading `X` matrix
                                          skip_assays = opt$corrected_only)

# check for corrected data 
corrected_assay_name <- paste(integration_method, "corrected", sep = "_")
if(!corrected_assay_name %in% assayNames(integrated_sce)){
  stop("integrated SCE missing corrected gene expression data in assays.")
}

# function to check for SVD or latent embeddings
check_pseudo_pca <- function(integrated_sce, 
                             pseudo_pca_name){
  # check corrected data and SVD were converted
  if(!pseudo_pca_name %in% reducedDimNames(integrated_sce)){
    stop(glue::glue("integrated SCE missing {pseudo_pca_name} in reduced dimensions slot."))
  }
}

# function to run UMAP with respective PCA name 
add_umap <- function(integrated_sce,
                     pseudo_pca_name,
                     umap_name){
  
  integrated_sce <- integrated_sce %>%
    scater::runUMAP(dimred = pseudo_pca_name,
                    name = umap_name)
}

# calculate UMAP dependent on integration method 
if(integration_method == "scanorama"){
  
  pseudo_pca_name <- "scanorama_SVD"
  check_pseudo_pca(integrated_sce,
                   pseudo_pca_name)
  
  umap_name <- "scanorama_UMAP"
  
  # run UMAP using the SVD as input 
  integrated_sce <- add_umap(integrated_sce,
                             pseudo_pca_name,
                             umap_name)
  
}

if(integration_method == "scvi"){
  
  pseudo_pca_name <- "scvi_latent"
  check_pseudo_pca(integrated_sce,
                   pseudo_pca_name)
  
  umap_name <- "scvi_UMAP"
  
  # run UMAP using the latent embeddings as input 
  integrated_sce <- add_umap(integrated_sce,
                             pseudo_pca_name,
                             umap_name)
  
}

# Write RDS --------------------------------------------------------------------

readr::write_rds(integrated_sce, opt$output_sce_file)
