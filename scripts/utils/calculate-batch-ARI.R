suppressPackageStartupMessages({
  library(magrittr)
  library(SingleCellExperiment)
})

#' Function to calculate batch ARI scores from an integrated SCE object
#'
#' This function uses a similar approach as used in this paper:
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
#' Data is subsampled to 80% of cells. K-means clustering is performed on a given
#' number of PCs (default 20), along a range of k's. BatchARI is then calculated.
#' This procedure is repeated 20x.
#'
#' @param integrated_sce The integrated SCE object
#' @param num_pcs The number of PCs to use during k-means clustering. Default: 20
#' @param seed Seed for initializing random sampling and k-means clustering
#' @param k_range Range of k values to use for k-means clustering. Default: `seq(5, 25, 5)`
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param unintegrated Indicates whether the provided data is intregated (`FALSE`; default) or
#'   integrated (`TRUE`).
#'   
#' @return Tibble with three columns: `rep`, representing the given downsampling replicate;
#'   `k`, the given k for k-means; `batch_ari`, the calculated batch ARI for the given
#'   combination of `rep` and `k`
calculate_batch_ari <- function(integrated_sce,
                                num_pcs = 20,
                                seed = NULL,
                                k_range = seq(5, 25, 5),
                                integration_method = NULL, 
                                unintegrated = FALSE) {

  # Set the seed for subsampling (_not bootstrap_; without replacement) and k-means
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
  all_pcs <- reducedDim(integrated_sce, reduced_dim_name)
  
  # Set up parameters
  frac_cells <- 0.8        # fraction of cells to downsample to
  nreps <- 20              # number of times to repeat sub-sampling procedure

  # perform calculations
  all_ari <- c()
  for (i in 1:nreps) {

    # Downsample PCs
    downsampled <- downsample_pcs_for_metrics(all_pcs, frac_cells, num_pcs)

    for (k in k_range) {

      # Perform k-means clustering
      clustering_result <- kmeans(x = downsampled$pcs, centers=k)
      downsampled_integrated_clusters <- unname(clustering_result$cluster)

      # Calculate batch ARI
      all_ari <- c(all_ari,
                   pdfCluster::adj.rand.index(downsampled_integrated_clusters, downsampled$batch_labels)
      )
    }
  }

  # Store results in a tibble
  ari_tibble <- tibble::tibble(
    rep = rep(1:nreps, each=length(k_range)),
    k   = rep(k_range, nreps),
    batch_ari = all_ari, 
    integration_method = integration_method
  )

  # Return tibble with batch ARI results which can further be summarized downstream
  return(ari_tibble)

}
