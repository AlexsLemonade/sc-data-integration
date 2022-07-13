library(SingleCellExperiment)

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
#' @param groups_to_integrate Array containing the covariates to consider during integration.
#' @param ... Additional parameters that may be passed to `harmony::HarmonyMatrix()`
#'
#' @return An integrated SCE object with two additional reducedDim fields: 
#'   `harmony_pcs` and `harmony_gene_matrix`
#' @examples 
#' integrate_harmony(combined_sce, c("batch))
integrate_harmony <- function(combined_sce, groups_to_integrate = c(), ...) {
  
  # Ensure columns given in groups_to_integrate are present in combined_sce ----
  if (length(groups_to_integrate) == 0) {
    stop("You must provide an array of groupings to integrate.")
  }
  if (!(all(groups_to_integrate %in% names(colData(combined_sce))))) {
    stop("The combined_sce object must contain covariate columns in colData.")
  }
  
  # Create metadata information for input to harmony ---------------------------
  harmony_metadata <- tibble::as_tibble(colData(combined_sce), 
                                        rownames = "cell_id") %>%
    dplyr::select(cell_id,
                  dplyr::all_of(groups_to_integrate))
  
  # Perform integration starting from PCs --------------------------------------
  # In this case, harmony uses pre-computed PCs during integration
  harmony_from_pcs <- harmony::HarmonyMatrix(
    data_mat  = reducedDim(combined_sce, "PCA"), 
    meta_data = harmony_metadata, 
    vars_use  = groups_to_integrate, 
    do_pca = FALSE, # We are passing in PCs
    ...
  )
  
  # Perform integration starting from the expression matrix --------------------
  # In this case, harmony computes PCs to use during integration
  harmony_from_gene_matrix <- harmony::HarmonyMatrix(
    data_mat  = logcounts(combined_sce), # Gene expression matrix
    meta_data = harmony_metadata, 
    vars_use  = groups_to_integrate, 
    do_pca = TRUE, # We are NOT passing in PCs
    ...
  )
  
  # Add new PCs from both approaches back into the combined_sce ----------------
  reducedDim(combined_sce, "harmony_pcs") <- harmony_from_pcs
  reducedDim(combined_sce, "harmony_gene_matrix") <- harmony_from_gene_matrix
  
  # Return the integrated SCE --------------------------------------------------
  return(combined_sce)
  
}



