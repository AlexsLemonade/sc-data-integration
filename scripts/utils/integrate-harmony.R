library(SingleCellExperiment)
library(magrittr) # pipe

#' Integrate a combined SCE with `harmony`
#' 
#' This function performs integration on a merged SCE object with different batches
#'  using the `harmony` R package. It performs integration using two approaches within
#'  `harmony`: i) integrating using provided PCs (here, those calculated by the 
#'  `scpca-downstream-analyses` workflow), and ii) integrating starting from the
#'  normalized gene matrix, meaning Harmony will calculate PCs to use for integration. 
#'  
#' @param combined_sce A combined SCE object. Must contain a cell column `batch`
#'   that indicates the different groups to be integrated.
#' @param groups_to_integrate Vector containing the covariates to consider during integration.
#' @param from_pca A logical indicating whether to integrate directly from PCs. Default: TRUE.
#' @param ... Additional parameters that may be passed to `harmony::HarmonyMatrix()`
#'
#' @return An integrated SCE object with the additional reducedDim field `harmony`
#'   representing the integrated PCs 
#' @examples 
#' integrate_harmony(combined_sce, c("batch"))
#' integrate_harmony(combined_sce, c("batch"), from_pca = FALSE) # start from gene expression matrix
integrate_harmony <- function(combined_sce, 
                              groups_to_integrate = c(), 
                              from_pca = TRUE,
                              ...) {
  
  # Perform checks ----------------------

  # Ensure groupings were provided
  if (length(groups_to_integrate) == 0) {
    stop("You must provide a vector of groupings to integrate.")
  }
  # Ensure groupings are present in the data
  if (!(all(groups_to_integrate %in% names(colData(combined_sce))))) {
    stop("The combined_sce object must contain covariate columns in colData.")
  }
  # Ensure PCs are present in the combined_sce object
  if (!("PCA" %in% reducedDimNames(combined_sce))) {
    stop("The combined_sce object must contain PCs.")
  }
  
  # Create metadata information for input to harmony ---------------------------
  harmony_metadata <- tibble::as_tibble(colData(combined_sce), 
                                        rownames = "cell_id") %>%
    dplyr::select(cell_id,
                  dplyr::all_of(groups_to_integrate))
  
  # Perform integration --------------------------------------
  if (from_pca) {
    # Perform integration using pre-computed PCs, as recommended
    harmony_results <- harmony::HarmonyMatrix(
      data_mat  = reducedDim(combined_sce, "PCA"), 
      meta_data = harmony_metadata, 
      vars_use  = groups_to_integrate, 
      do_pca = FALSE, # We are passing in PCs
      ...
    )
  } else {
    # Perform integration starting from the expression matrix
    harmony_results<- harmony::HarmonyMatrix(
      data_mat  = logcounts(combined_sce), # Gene expression matrix
      meta_data = harmony_metadata, 
      vars_use  = groups_to_integrate, 
      do_pca = TRUE, # We are NOT passing in PCs
      ...
    )
  }
  
  # Add new PCs back into the combined_sce ----------------
  reducedDim(combined_sce, "harmony") <- harmony_results

  # Return the integrated SCE --------------------------------------------------
  return(combined_sce)
  
}



