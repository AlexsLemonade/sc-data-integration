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
#' @param num_pcs Number of PCs to consider during regression calculation. Default: 50
#' @param significance_threshold Threshold for considering a corrected P-value from PC~batch
#'  regression as significant. Default: 0.05
#' @param integration_method The name of the method that was used for integration
#'  to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param unintegrated Indicates whether the provided data is intregated (`FALSE`; default) or
#'   integrated (`TRUE`).
#' @param seed Seed for initializing random sampling
#'
#' @return Tibble with four columns: `rep`, the given downsampling replicate; 
#' `integration_method`, the given integration method (`NA` if unintegrated); 
#' `batch_variance`, the  total contribution of the batch effect to the variance; 
#' `pc_regression_scaled`, proxy for batch effect
#'
calculate_pca_regression <- function(integrated_sce,
                                     batch_column = "batch",
                                     num_pcs = 50,
                                     significance_threshold = 0.05,
                                     integration_method = NULL, 
                                     unintegrated = FALSE,
                                     seed = NULL) {

  # Set seed if provided
  set.seed(seed)

  # Settings depending on whether data is integrated or not
  if (unintegrated){
    # In the end, we'll return NA in the data frame for integration method
    integration_method <- NA
    
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

  
  # Set up parameters - using the same as for ARI and AWS
  frac_cells <- 0.8        # fraction of cells to downsample to
  nreps <- 20              # number of times to repeat sub-sampling procedure

  batch_variances <- c()
  pc_reg_scales   <- c()
  
  # Perform over `nreps` subsets of data 
  for (i in 1:nreps) {

    # Prepare data for modeling, beginning with downsampling the PCs (keep all PCs, but only frac of cells)
    pc_tibble <- downsample_pcs_for_metrics(pcs, frac_cells, ncol(pcs))$pcs %>%
      tibble::as_tibble(pcs, rownames = "cellname") %>%
      # Pull out batch column to regress against
      tidyr::separate(cellname, into = c("barcode", "batch"), sep = "-") %>%
      dplyr::select(-barcode) %>%
      # Create long tibble with columns: batch, PC_index (<1:num_pcs>), and PC (the PCs)
      tidyr::pivot_longer(-batch, names_to = "PC_index", values_to = "PC") %>%
      dplyr::mutate(PC_index = as.numeric(
        stringr::str_replace(PC_index, "^V", ""))
      )
    
    
    # For each PC, regress against `batch` and obtain R^2 and associated P-value 
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
    explained_variance <- (pcs_variance / sum_of_variance)
  
    # This is the approximation of: "the total contribution of the batch effect to the variance in the data" (pg 50, Methods section)
    batch_variance <- sum(pc_regression_results$r.squared*explained_variance)
  
  
    # Calculate the "sum of explained variance of all PCs with significant <R^2 from the batch regression>
    #  scaled by the variance explained by the top `num_pcs` PCs as a proxy for the batch effect" (pg 50, Methods section)
    
    # First, determine which corrected P-values are significant
    corrected_pvalues_sig <- unname(p.adjust(pc_regression_results$p.value, method = 'BH')) < significance_threshold
  
    # keep only top `num_pcs`
    corrected_pvalues_sig <- corrected_pvalues_sig[1:num_pcs]
    
    # This the "proxy for the batch effect, scaled 0-1 where 0 = no effect and 1 = strong effect
    # Perform across top `num_pcs` only
    pc_reg_scale <- sum(explained_variance[corrected_pvalues_sig])/sum(explained_variance[1:num_pcs])
  
    
    # Store quantities
    batch_variances <- c(batch_variances, batch_variance)
    pc_reg_scales   <- c(pc_reg_scales, pc_reg_scale)
  
  }
  
  
  # Create and return result tibble -----
  results <- tibble::tibble(
    rep = 1:nreps,
    integration_method = integration_method,
    pc_batch_variance = batch_variances,
    pc_regression_scaled = pc_reg_scales,
    .name_repair = 'minimal'
  )
  
  return(results)
    
}
