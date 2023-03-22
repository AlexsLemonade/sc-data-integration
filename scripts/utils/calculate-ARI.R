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
  
  # Pull out the PCs or analogous reduction, and remove batch NAs along the way
  pcs <- reducedDim(integrated_sce, reduced_dim_name)
  # we do not need the indices here, so just save the pcs directly
  final_pcs <- remove_batch_nas_from_pcs(pcs, colData(integrated_sce)[,batch_column])[["pcs"]]

  # Set up parameters
  frac_cells <- 0.8        # fraction of cells to downsample to
  nreps <- 20              # number of times to repeat sub-sampling procedure

  # perform calculations
  all_ari <- c()
  for (i in 1:nreps) {

    # Downsample PCs
    downsampled <- downsample_pcs_for_metrics(final_pcs, frac_cells, num_pcs)

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

#' Calculate Reverse ARI to compare clustering pre and post integration
#'
#' @param individual_sce_list List of SCE objects containing all batches to use for 
#'   calculating reverse ARI. Must be the original, pre-merged object.
#'   The list must be a named list with the batch names.
#' @param integrated_sce The integrated SCE object
#' @param seed Seed for initializing random sampling and k-means clustering
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param unintegrated Indicates whether the provided data is integrated (`FALSE`; default) or
#'   integrated (`TRUE`).
#' @param batch_column The variable in `integrated_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "batch".
#'
#' @return Data frame with three columns: `ari`, the calculated ARI, 
#'   the corresponding `batch_id` found in the `batch_column` 
#'   for the original SCE object, and the `integration_method`, the provided integration method
#'
calculate_reverse_ari <- function(individual_sce_list,
                                  integrated_sce,
                                  seed = NULL,
                                  integration_method = NULL,
                                  unintegrated = FALSE,
                                  batch_column = "batch"){
  set.seed(seed)
  
  # check that list of SCE objects is named
  batch_ids <- names(individual_sce_list)
  if(is.null(batch_ids)){
    stop("Must provide a named list of SCE objects.")
  } else {
    # make sure the batch ids provided match between the list and the integrated object 
    if(!all(batch_ids %in% unique(colData(integrated_sce)[, batch_column]))){
      stop("Names of provided SCE objects included in the individual SCE object 
           do not match batch IDs present in the batch_column of the integrated object")
    }
  }
  
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

  
  # get pcs for individual objects
  ind_pcs <- purrr::map(individual_sce_list, get_individual_pcs)
  
  # Pull out the PCs or analogous reduction, and remove batch NAs along the way
  integrated_pcs <- reducedDim(integrated_sce, reduced_dim_name)
  
  # cluster integrated pcs only one time
  integrated_clustering_result <- bluster::clusterRows(integrated_pcs, 
                                                       bluster::NNGraphParam(cluster.fun = "louvain",
                                                                             type = "jaccard",
                                                                             k = k)) |> 
    set_names(rownames(integrated_pcs)) # make sure to set the names with the batch ids
  
  
  # for every batch id, cluster and then calculate ari for that batch 
  all_ari <- batch_ids |> 
    purrr::map_chr(\(batch){
      
      ind_clustering_result <- bluster::clusterRows(ind_pcs[[batch]],
                                                    bluster::NNGraphParam(cluster.fun = "louvain",
                                                                          type = "jaccard",
                                                                          k = k))
      # extract clusters from integrated clustering for batch 
      clusters_to_keep <- grep(batch, names(integrated_clustering_result))
      batch_integrated_clusters <- integrated_clustering_result[clusters_to_keep]
      
      ari <- bluster::pairwiseRand(ind_clustering_result, 
                                   batch_integrated_clusters, 
                                   mode = "index")
      
    })
  
  # create data frame with ari and batch id
  ari_df <- data.frame(
    ari = all_ari,
    batch_id = batch_ids,
    integration_method = integration_method
  )

  return(ari_df)
}
