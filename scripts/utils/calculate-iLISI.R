suppressPackageStartupMessages({
  library(lisi)
  library(magrittr)
  library(SingleCellExperiment)
})

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
  all_integration_methods <- c("fastMNN", "harmony", "rpca", "cca", "scvi", "scanorama")
  if (is.null(integration_method)) {
    stop("An `integration_method` must be provided.")
  } else {
    integration_method <- tolower(integration_method)
    if (!(integration_method %in% tolower(all_integration_methods))) {
      stop(
        paste("The `integration_method` must be one of: ",
              paste(all_integration_methods, collapse = ", ")
        )
      )
    }
  }
  
  # Create data frame with batch information to provide to `compute_lisi()`
  batch_df <- data.frame(batch = colData(integrated_sce)[,batch_column])
  
  # Extract the PCs to provide to `compute_lisi()`
  pcs <- reducedDim(integrated_sce, 
                    paste0(integration_method, "_PCA"))
  

  
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
    dplyr::mutate(cell = colnames(integrated_sce),
                  integration_method = integration_method) %>%
    # split cell into cell_barcode and library
    tidyr::separate(cell, 
                    into = c("cell_barcode", "library"), 
                    sep = "-")

  
  # Return the tibble
  return(lisi_result)

}