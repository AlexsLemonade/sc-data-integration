#' Integrate a combined SCE with Harmony
#' 
#' This function performs integration on a merged SCE object with different batches
#'  using the harmony R package. It performs integration using two approaches within
#'  Harmony: i) integrating using provided PCs (here, those calculated by the 
#'  `scpca-downstream-analyses` workflow), and ii) integrating starting from the
#'  normalized gene matrix, meaning Harmony will calculate PCs to use for integration. 
#'  
#' @param combined_sce A combined SCE object. Must contain a cell column `batch`
#'   that indicates the different groups to be integrated.
#'
#' @return An integrated SCE object with two additional reducedDim names: 
#'   `harmony_pcs` and `harmony_gene_matrix`
integrate_harmony <- function(combined_sce) {
  
  # Ensure `batch` column is present in combined_sce ---------------------------
  if (!("batch" %in% names(colData(combined_sce)))) {
    stop("The combined_sce object must contain a cell data field `batch` to perform integration.")
  }
  
  # Create metadata information for input to Harmony ---------------------------
  harmony_metadata <- tibble::tibble(
    cell_id = rownames(colData(combined_sce)),
    batch = colData(combined_sce)$batch
  ) 
  
  # Perform integration starting from PCs --------------------------------------
  # In this case, Harmony uses pre-computed PCs during integration
  harmony_from_pcs <- HarmonyMatrix(
    data_mat  = reducedDim(combined_sce, "PCA"), # Matrix with coordinates for each cell (row) along many PCs (columns)
    meta_data = harmony_metadata, 
    vars_use  = "batch", # column in meta_data that indicates groups to integrate
    do_pca = FALSE # We are passing in PCs
  )
  
  # Perform integration starting from the expression matrix --------------------
  # In this case, Harmony computes PCs to use during integration
  harmony_from_gene_matrix <- HarmonyMatrix(
    data_mat  = logcounts(combined_sce), # Gene expression matrix
    meta_data = harmony_metadata, # Dataframe with information for each cell (row)
    vars_use  = "batch",
    do_pca = TRUE # We are NOT passing in PCs
  )
  
  # Add new PCs from both approaches back into the combined_sce ----------------
  reducedDim(combined_sce, "harmony_pcs") <- harmony_from_pcs
  reducedDim(combined_sce, "harmony_gene_matrix") <- harmony_from_gene_matrix
  
  # Return the integrated SCE --------------------------------------------------
  return(combined_sce)
  
}



