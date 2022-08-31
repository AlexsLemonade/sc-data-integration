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
#' @param unintegrated Indicates whether the provided data is integrated (`FALSE`; default) or
#'   integrated (`TRUE`).
#'   
#' @return Tibble with five columns with one row per cell. Columns are `ilisi_score`, 
#'   `cell_name` (combined barcode and library), `cell_barcode`, `library` and 
#'   `integration_method`
calculate_ilisi <- function(integrated_sce,
                            batch_column = "batch", 
                            integration_method = NULL, 
                            unintegrated = FALSE) {
  
  # Settings depending on whether data is integrated or not
  if (unintegrated){
    # In the end, we'll return "unintegrated" in the data frame for integration method
    integration_method <- "unintegrated"
    
    # use simply "PCA" for reduced dimensions
    reduced_dim_name <- "PCA"
  } else {
    # Check integration method
    integration_method <- check_integration_method(integration_method)
    
    # Get name for reduced dimensions
    reduced_dim_name <- get_reduced_dim_name(integration_method)
  }

  # Pull out the PCs or analogous reduction
  pcs <- reducedDim(integrated_sce, reduced_dim_name)
  
  # Create data frame with batch information to provide to `compute_lisi()`
  batch_df <- data.frame(batch = colData(integrated_sce)[,batch_column])
  
  
  # `lisi_result` is a tibble with per-cell scores, the score roughly means:
  #   "how many different categories are represented in the local neighborhood of the given cell?"
  #   With 2 batches, then, values close to 2 indicate good integration
  lisi_result <- lisi::compute_lisi(pcs, 
                                    # define the batches
                                    batch_df, 
                                    # which variables in `batch_df` to compute lisi for
                                    "batch") %>% 
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
