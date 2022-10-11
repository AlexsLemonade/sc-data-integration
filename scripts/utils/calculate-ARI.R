suppressPackageStartupMessages({
  library(magrittr)
  library(SingleCellExperiment)
})

#' Function to calculate ARI scores from an integrated SCE object
#'
#' This function uses a similar approach as used in this paper:
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
#' Data is subsampled to 80% of cells. K-means clustering is performed on a given
#' number of PCs (default 20), along a range of k's. ARI (either batch or celltype) 
#' is then calculated. This procedure is repeated 20x.
#'
#' @param integrated_sce The integrated SCE object
#' @param num_pcs The number of PCs to use during k-means clustering. Default: 20
#' @param seed Seed for initializing random sampling and k-means clustering
#' @param k_range Range of k values to use for k-means clustering. Default: `seq(5, 25, 5)`
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param unintegrated Indicates whether the provided data is integrated (`FALSE`; default) or
#'   integrated (`TRUE`).
#' @param batch_column The variable in `integrated_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "batch".
#'   
#' @return Tibble with four columns: `rep`, representing the given downsampling replicate;
#'   `k`, the given k for k-means; `ari`, the calculated ARI for the given
#'   combination of `rep` and `k`; the `integration_method`, the given integration method
calculate_ari <- function(integrated_sce,
                          num_pcs = 20,
                          seed = NULL,
                          k_range = seq(5, 25, 5),
                          integration_method = NULL, 
                          unintegrated = FALSE, 
                          batch_column = "batch") {
  # Set the seed for subsampling (_not bootstrap_; without replacement) and k-means
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
  }
  
  # Pull out the PCs or analogous reduction
  all_pcs <- reducedDim(integrated_sce, reduced_dim_name)

  
  # We need to _remove columns_ (cells) with unknown batch_column. 
  # This will usually mean removing NA cell types in the context of celltype ARI calculations
  batches <- colData(integrated_sce)[,batch_column]
  retain_indices <- which(!is.na(batches))
  batches <- batches[retain_indices]
  all_pcs <- all_pcs[retain_indices,]
  
  # check dimensions still match:
  if (nrow(all_pcs) != length(batches)) {
    stop("Incompatable PC and batch information dimensions after removing NAs.")
  }
  
  # Set PC rownames to be the batches for obtaining downsampled labels
  rownames(all_pcs) <- batches
  
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
    ari = all_ari, 
    integration_method = integration_method
  )

  # Return tibble with batch ARI results which can further be summarized downstream
  return(ari_tibble)

}
