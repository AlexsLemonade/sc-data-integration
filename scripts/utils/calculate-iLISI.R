suppressPackageStartupMessages({
  library(lisi)
  library(magrittr)
  library(SingleCellExperiment)
})

# script with `check_integration_method()`
source(
  here::here(
    "scripts", 
    "utils", 
    "integration-helpers.R"
  )
)

#' Function to calculate iLISI scores from an integrated SCE object
#'
#' @param integrated_sce The integrated SCE object
#' @param batch_column The variable in `integrated_sce` indicating batches. Default
#'   is "batch".
#' @param integration_method The name of the method that was used for integration 
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#'
#' @return Tibble with four columns with one row per cell. Columns are `ilisi_score`, 
#'   `cell_barcode`, `library` and `integration_method`
calculate_ilisi <- function(integrated_sce,
                            batch_column = "batch", 
                            integration_method = NULL) {
  
  # Check integration method
  integration_method <- check_integration_method(integration_method)
  
  # Create data frame with batch information to provide to `compute_lisi()`
  batch_df <- data.frame(batch = colData(integrated_sce)[,batch_column])
  
  
  # Extract the PCs (or similar) to provide to `compute_lisi()`
  
  # If we are using either scanorama or scvi, there will be a different name:
  if (integration_method == "scanorama") {
    reduced_dim_name <- "scanorama_SVD"
  } else if (integration_method == "scvi") {
    reduced_dim_name <- "scvi_latent"
  } else {
    reduced_dim_name <- paste0(integration_method, "_PCA")
  }
  pcs <- reducedDim(integrated_sce, reduced_dim_name)
  
  # `lisi_result` is a tibble with per-cell scores, the score roughly means:
  #   "how many different categories are represented in the local neighborhood of the given cell?"
  #   With 2 batches, then, values close to 2 indicate good integration
  lisi_result <- lisi::compute_lisi(pcs, 
                                    # define the batches
                                    batch_df, 
                                    # which variables in `batch_df` to compute lisi for
                                    batch_column) %>% 
    tibble::as_tibble() %>%
    # Rename the result column to `ilisi_score`
    dplyr::rename(ilisi_score = batch) %>%
    # Add in the cell, library ID, and integration method
    dplyr::mutate(cell_name = colnames(integrated_sce),
                  integration_method = integration_method) %>%
    # split cell into cell_barcode and library
    tidyr::separate(cell_name, 
                    into = c("cell_barcode", "library"), 
                    sep = "-", 
                    # keep the original cell!
                    remove = FALSE) 

  
  # Return the tibble
  return(lisi_result)

}
