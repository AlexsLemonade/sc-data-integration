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
#' @param num_pcs Number of PCs to consider during kBET calculation. Default: 20
#' @param k0_fraction_range Range of fractions of the sample size to set k0 (neighborhood
#'  size) to when running `kBET`. Default values are based on the following paper:
#'    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
#      "Following the example in the kBET paper, we chose the k input value equal to
#       5%, 10%, 15%, 20%, and 25% of the sample size"
#' @param unintegrated Indicates whether the provided data is integrated (`FALSE`; default) or
#'   integrated (`TRUE`).
#' @param seed Random seed to set for kBET calculations
#'
#' @return Tibble with five columns: `k0_fraction`, the fraction of total sample sized used to
#'  assign the neighborhood size used by `kBET`; `rep`, the given subsample for that `k0`;
#'  `integration_method`,the given integration method; `kbet_stat_type`, either 
#'  "observed_rejection_rate" or "expected_rejection rate"; `kbet_stat`, the given rejection
#'  rate value seen in `kbet_stat_type` for the given subsample/k0 combination
#'  
calculate_kbet <- function(integrated_sce,
                           batch_column = "batch",
                           integration_method = NULL,
                           unintegrated = FALSE,
                           num_pcs = 20,
                           k0_fraction_range = c(0.05, 0.10, 0.15, 0.20, 0.25),
                           seed = NULL) {

  # Set seed if given
  set.seed(seed)


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
    
    # if integration method is scvi, make sure num pcs is set to 10 or less
    if(integration_method == "scvi"){
      if(num_pcs > 10){
        num_pcs = 10
      }
    }
  }
  
  # Pull out the PCs or analogous reduction
  pcs <- reducedDim(integrated_sce, reduced_dim_name)
  
  # Calculate the mean batch size for k0 calculation
  batches <- colData(integrated_sce)[,batch_column]
  mean_batch_size <- table(batches) %>%
    tibble::as_tibble() %>%
    dplyr::summarize(mean_batch_size = mean(n)) %>%
    dplyr::pull(mean_batch_size)

  # Run kBET across given range of k0
  k0_range <- ceiling(mean_batch_size * k0_fraction_range)

  kbet_observed_rates <- c()
  kbet_expected_rates <- c()
  for (k0 in k0_range) {
    kbet_result <- kBET(pcs,
                        batches,
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

      kbet_observed_rates <- c(kbet_observed_rates,
                                # observed (actual) rejection rates for each of the 100 subsamples
                                kbet_result$stats$kBET.observed)
      kbet_expected_rates <- c(kbet_expected_rates,
                               # expected rejection rates for each of the 100 subsamples
                               kbet_result$stats$kBET.expected)

  }

  # Place results into tibble
  kbet_results <- tibble::tibble(
    k0_fraction = rep(k0_fraction_range, each = 100),
    rep = rep(1:100, length(k0_range)),
    observed_rejection_rate = kbet_observed_rates,
    expected_rejection_rate = kbet_expected_rates,
    integration_method = integration_method
  ) %>%
    # tidy it up a bit for easier plotting
    tidyr::pivot_longer(
      dplyr::ends_with("rejection_rate"), 
      names_to = "kbet_stat_type",
      values_to = "kbet_stat"
    )
    


  # Return results
  return(kbet_results)

}
