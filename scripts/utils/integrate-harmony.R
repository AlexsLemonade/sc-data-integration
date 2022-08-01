library(SingleCellExperiment)
library(magrittr) 

source(
  here::here(
    "scripts", 
    "utils",
    "integration-helpers.R"
  )
)


#' Integrate a combined SCE with `harmony`
#' 
#' This function performs integration on a merged SCE object with different batches
#'  using the `harmony` R package. It can perform integration using one of two approaches within
#'  `harmony`: i) integrating using provided PCs (here, those calculated by the 
#'  `scpca-downstream-analyses` workflow, this is the default approach), and ii) integrating starting from the
#'  normalized gene matrix, meaning Harmony will calculate PCs to use for integration. 
#'  
#' @param combined_sce A combined SCE object. Must contain a cell column `batch`
#'   that indicates the different groups to be integrated.
#' @param batch_column The variable in `combined_sce` indicating batches. Default 
#'   is "batch".
#' @param covariate_cols Optional vector containing additional covariates to consider during integration.
#' @param from_pca A boolean indicating whether to integrate directly from PCs. Default: TRUE.
#' @param seed Random seed to set for `harmony` integration. A seed will only
#'  be set if this is not `NULL` (the default).
#' @param return_uncorrected_expression Boolean indicating whether to return the
#'  uncorrected expression values in the returned SCE object. Default: `FALSE`.
#' @param ... Additional parameters that may be passed to `harmony::HarmonyMatrix()`
#'
#' @return An integrated SCE object with the additional reducedDim field `harmony`
#'   representing the integrated PCs 
#' @examples 
#' integrate_harmony(combined_sce)
#' integrate_harmony(combined_sce, covariate_cols = "sample", from_pca = FALSE) # start from gene expression matrix
integrate_harmony <- function(combined_sce, 
                              batch_column = "batch",
                              covariate_cols = c(), 
                              from_pca = TRUE,
                              seed = NULL,
                              return_uncorrected_expression = FALSE,
                              ...) {
  
  if (!(is.null(seed))) {
    set.seed(seed)
  }
  
  # Perform checks ----------------------

  # Ensure batch column is present
  if (!(batch_column %in% names(colData(combined_sce)))) {
    stop("The provided `batch_column` column must be in `combined_sce` colData.")
  }
  
  # Ensure covariates are present in the data, if provided
  if (length(covariate_cols) > 0) {
    if (!(all(covariate_cols %in% names(colData(combined_sce))))) {
      stop("The combined_sce object must contain the given covariate columns in colData.")
    }
  }
  
  # Ensure PCs are present in the combined_sce object
  if (from_pca && !("PCA" %in% reducedDimNames(combined_sce))) {
    stop("The combined_sce object must contain PCs.")
  }
  
  # Create metadata information for input to harmony ---------------------------
  # First, combine together all the covariates with batch to use together:
  covariate_cols <- c(batch_column, covariate_cols) 
  harmony_metadata <- tibble::as_tibble(colData(combined_sce), 
                                        rownames = "cell_id") %>%
    dplyr::select(cell_id,
                  dplyr::all_of(covariate_cols))
  
  # Perform integration --------------------------------------
  if (from_pca) {
    # Perform integration using pre-computed PCs, as recommended
    harmony_results <- harmony::HarmonyMatrix(
      data_mat  = reducedDim(combined_sce, "PCA"), 
      meta_data = harmony_metadata, 
      vars_use  = covariate_cols, 
      do_pca = FALSE, # We are passing in PCs
      ...
    )
  } else {
    # Perform integration starting from the expression matrix
    harmony_results<- harmony::HarmonyMatrix(
      data_mat  = logcounts(combined_sce), # Gene expression matrix
      meta_data = harmony_metadata, 
      vars_use  = covariate_cols, 
      do_pca = TRUE, # We are NOT passing in PCs
      ...
    )
  }
  
  # Add new PCs back into the combined_sce ----------------
  reducedDim(combined_sce, "harmony_PCA") <- harmony_results
  
  # Perform UMAP with the new PCs --------------------------
  combined_sce <- perform_dim_reduction(combined_sce, prefix = "harmony")
  
  # Remove uncorrected expression values, unless otherwise specified ----
  if (!return_uncorrected_expression) {
    combined_sce <- remove_uncorrected_expression(combined_sce)
  }
  
  # Return the integrated SCE --------------------------------------------------
  return(combined_sce)
  
}



