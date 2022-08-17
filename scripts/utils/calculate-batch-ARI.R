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
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#'
#' @return Tibble with three columns: `rep`, representing the given downsampling replicate;
#'   `k`, the given k for k-means; `batch_ari`, the calculated batch ARI for the given
#'   combination of `rep` and `k`
calculate_batch_ari <- function(integrated_sce,
                                num_pcs = 20,
                                seed = NULL,
                                integration_method = NULL) {

  # Check integration method
  integration_method <- check_integration_method(integration_method)

  # Set the seed for subsampling (_not bootstrap_; without replacement) and k-means
  set.seed(seed)

  # Set up parameters
  frac_cells <- 0.8        # fraction of cells to downsample to
  k_range <- seq(5, 25, 5) # range of k for k-means clustering
  nreps <- 20              # number of times to repeat sub-sampling procedure

  # pull out relevant information for calculations
  all_integrated_pcs <- reducedDim(integrated_sce, paste0(integration_method, "_PCA"))
  num_cells <- nrow(all_integrated_pcs)


  # perform calculations
  all_ari <- c()
  for (i in 1:nreps) {

    # Get downsampling (without replacement) indices
    downsampled_indices <- sample(1:num_cells,
                                  frac_cells*num_cells,
                                  replace = FALSE)

    # Extract PCs for downsample, considering only the top `num_pcs`
    downsampled_integrated_pcs <- all_integrated_pcs[downsampled_indices,1:num_pcs]

    for (k in k_range) {

      # Perform k-means clustering
      clustering_result <- kmeans(x = downsampled_integrated_pcs, centers=k)
      downsampled_integrated_clusters <- unname(clustering_result$cluster)

      # Obtain batch labels as integer values
      downsampled_batch_labels <- stringr::str_replace_all(
        rownames(downsampled_integrated_pcs),
        "^[ACGT]+-",
        ""
      ) %>%
        as.factor() %>%
        as.numeric()


      # Calculate batch ARI
      all_ari <- c(all_ari,
                   pdfCluster::adj.rand.index(downsampled_integrated_clusters, downsampled_batch_labels)
      )
    }
  }

  # Store results in a tibble
  ari_tibble <- tibble::tibble(
    rep = rep(1:nreps, each=length(k_range)),
    k   = rep(k_range, nreps),
    batch_ari = all_ari
  )

  # Return tibble with batch ARI results which can further be summarized downstream
  return(ari_tibble)

}
