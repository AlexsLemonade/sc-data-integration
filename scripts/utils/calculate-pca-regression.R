suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

source(
  here::here(
    "scripts",
    "utils",
    "integration-helpers.R"
  )
)

#' Function to perform PCA regression for batch effect evaluation as described 
#'   in the `kBET` manuscript
#'
#' @param integrated_sce The integrated SCE object
#' @param batch_column The variable in `integrated_sce` indicating batches. Default
#'   is "batch".
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param num_pcs Number of PCs to consider during regression calculation. Default: 50
#' @significance_threshold Threshold for considering a corrected P-value from PC~batch 
#'  regression as significant. Default: 0.05

#'
#' @return Tibble with five columns: `k0_fraction`, the fraction of total sample sized used to
#'  assign the neighborhood size used by `kBET`; `rep`, the given subsample for that `k0`;
#'  `integration_method`,the given integration method; `kbet_stat_type`, either 
#'  "observed_rejection_rate" or "expected_rejection rate"; `kbet_stat`, the given rejection
#'  rate value seen in `kbet_stat_type` for the given subsample/k0 combination
#'  

calculate_pca_regression <- function(integrated_sce, 
                                     batch_column = "batch", 
                                     num_pcs = 50,
                                     significance_threshold = 0.05,
                                     integration_method = NULL) {
  
  
  # Check integration method
  integration_method <- check_integration_method(integration_method)
  
  # Get PCs or analagous reduced dimensions
  reduced_dim_name <- get_reduced_dim_name(integration_method)
  pcs <- reducedDim(integrated_sce, reduced_dim_name)
  
  # Subset to given num_pcs PCs
  pca.data <- pcs[,1:num_pcs]
  
  
  perform_regression <- function(df) {
    # Adapted from https://github.com/theislab/kBET/blob/f35171dfb04c7951b8a09ac778faf7424c4b6bc0/R/kBET-utils.R#L168-L181
    # Note that we do not return this quantity: https://github.com/theislab/kBET/blob/f35171dfb04c7951b8a09ac778faf7424c4b6bc0/R/kBET-utils.R#L177
    #   because it appears to be a hold-over from this line: https://github.com/theislab/kBET/blob/f35171dfb04c7951b8a09ac778faf7424c4b6bc0/R/kBET-utils.R#L160
    #   , and it doesn't make sense to retain only one coefficient's of the P-values; the P-value we want is the REGRESSION P-value, as described in the paper.
    #   It also appears the code in their main regression function carries over a bug by using that P-value as the regression P-value, when it is not:
    #    https://github.com/theislab/kBET/blob/master/R/pcRegression.R
    #   The paper further describes usage of an F-statistic. While their function does return the F stat's P-value, that P-value is never used in their code.
    #    Therefore we don't calculate it here.
    reg_result <- lm(PC ~ batch, data = df) %>%
      broom::glance() %>%
      dplyr::select(r.squared, p.value)
    return(reg_result)
  }
  
  # For each PC, regress against `batch` and obtain R^2 and associated P-value ----------
  
  # Prepare data for modeling
  pc_tibble <- tibble::as_tibble(pca.data, rownames = "cellname") %>%
    # Pull out batch column to regress against
    tidyr::separate(cellname, into = c("barcode", "batch"), sep = "-") %>%
    dplyr::select(-barcode) %>%
    # Create long tibble with columns: batch, PC_index (<1:num_pcs>), and PC (the PCs)
    tidyr::pivot_longer(-batch, names_to = "PC_index", values_to = "PC") %>%
    dplyr::mutate(PC_index = as.numeric(
      stringr::str_replace(PC_index, "^V", ""))
    )
  
  # Perform regressions
  pc_regression_results <- pc_tibble %>%
    # Create nested data for each PC
    dplyr::group_by(PC_index) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    # Perform regression PC~batch and unnest
    dplyr::mutate(fit = purrr::map(data, perform_regression)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(fit) 
  
  
  # Calculate metrics from regression  results --------------
  # All calculations are adapted from https://github.com/theislab/kBET/blob/master/R/pcRegression.R
  
  # Calculate the "total contribution of the batch effect to the variance"
  pcs_variance <- pc_tibble %>%
    dplyr::group_by(PC_index) %>% 
    dplyr::summarize(pc_variance = sd(PC)**2) %>%
    dplyr::pull(pc_variance)
  
  # Calculate \sum variance and the percentage explained by PCs 
  sum_of_variance <- sum(pcs_variance) 
  explained_variance <- (pcs_variance / sum_of_variance)/100
   
  # This is the approximation of: "the total contribution of the batch effect to the variance in the data" (pg 50, Methods section)
  batch_variance <- sum(pc_regression_results$r.squared*explained_variance) 
  
  
  # Calculate the "sum of explained variance of all PCs with significant <R^2 from the batch regression>
  #  scaled by the variance explained by the top `num_pcs` PCs as a proxy for the batch effect" (pg 50, Methods section)
  
  # First, determine which corrected P-values are significant  
  corrected_pvalues_sig <- unname(p.adjust(pc_regression_results$p.value, method = 'BH')) < significance_threshold
  
  # This the "proxy for the batch effect, scaled 0-1 where 0 = mo effect and 1 = strong effect
  pc_reg_scale <- sum(explained_variance[corrected_pvalues_sig])/sum(explained_variance) 
  
  
  
  # Return a 1-row tibble -----
  return(
    tibble::tibble(
      integration_method = integration_method,
      pc_batch_variance = batch_variance,
      pc_regression_scaled = pc_reg_scale
    )
  )
  
}