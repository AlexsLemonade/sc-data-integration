suppressPackageStartupMessages({
  library(kBET)
  library(SingleCellExperiment)
})

source(
  here::here(
    "scripts", 
    "utils", 
    "integration-helpers.R"
  )
)

#' Function to calculate kBET rejection rates from an integrated SCE object
#'
#' @param integrated_sce The integrated SCE object
#' @param batch_column The variable in `integrated_sce` indicating batches. Default
#'   is "batch".
#' @param integration_method The name of the method that was used for integration 
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param num_pcs Number of PCs to consider during kBET calculation. Default: 20 (TODO!)
#' @param k0_fraction_range Range of fractions of the sample size to set k0 (neighborhood 
#'  size) to when running `kBET`. Values are based on the following paper:
#'    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
#      "Following the example in the kBET paper, we chose the k input value equal to 
#       5%, 10%, 15%, 20%, and 25% of the sample size"
#' @param seed Random seed to set for kBET calculations
#'
#' @return Tibble with four columns: `k0_fraction`, the fraction of total sample sized used to 
#'  assign the neighborhood size used by `kBET`; `rep`, the given subsample for that `k0`; 
#'  `kbet_rejection_rate`, the `kBET` rejection rate for that subsample; 
#'  `integration_method`,the given integration method
#'  
calculate_kbet <- function(integrated_sce,
                           batch_column = "batch", 
                           integration_method = NULL,
                           num_pcs = 20,
                           k0_fraction_range = c(0.05, 0.10, 0.15, 0.20, 0.25),
                           seed = NULL) {
  
  # Set seed if given
  set.seed(seed)
  
  
  # Check integration method
  integration_method <- check_integration_method(integration_method)
  
  
  # Get PCs or analagous reduced dimensions
  reduced_dim_name <- get_reduced_dim_name(integration_method)
  pcs <- reducedDim(integrated_sce, reduced_dim_name)

  
    k0_range <- ceiling(nrow(pcs) * k0_fraction_range)
  kbet_rejection_rates <- c()
  for (k0 in k0_range) {
    kbet_result <- kBET(pcs,
                        colData(integrated_sce)[,batch_column], 
                        # we are providing PCs directly
                        do.pca = FALSE,
                        # specify number of PCs
                        dim.pca = num_pcs,
                        # Don't make a boxplot
                        plot = FALSE,
                        # use the provided neighborhood size
                        k0 = k0, 
                        # prevent kBET from changing the provided neighborhood size
                        heuristic = FALSE)
    
      kbet_rejection_rates <- c(kbet_rejection_rates, 
                                # rejection rates for each of the 100 subsamples
                                kbet_result$stats$kBET.observed)

  }
  
  # Place results into tibble
  kbet_results<- tibble::tibble(
    k0_fraction = rep(k0_fraction_range, each = 100),
    rep = rep(1:100, length(k0_range)),
    kbet_rejection_rate = kbet_rejection_rates,
    integration_method = integration_method
  )

  
  # Return results
  return(kbet_results)
  
}
