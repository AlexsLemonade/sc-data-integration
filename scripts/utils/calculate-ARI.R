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
#' @return Tibble with five columns: `rep`, representing the given downsampling replicate;
#'   `k`, the given k for k-means; `ari`, the calculated ARI for the given
#'   combination of `rep` and `k`; the corresponding `batch_id` found in the `batch_column` 
#'   for the original SCE object, and the `integration_method`, the provided integration method
#'
calculate_reverse_ari <- function(individual_sce_list,
                                  integrated_sce,
                                  num_pcs = 20,
                                  seed = NULL,
                                  k_range = seq(5, 25, 5),
                                  integration_method = NULL,
                                  unintegrated = FALSE,
                                  batch_column = "batch"){
  set.seed(seed)
  
  if(is.null(names(individual_sce_list))){
    stop("Must provide a named list of SCE objects.")
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
  
  # Set up parameters
  frac_cells <- 0.8        # fraction of cells to downsample to
  nreps <- 20              # number of times to repeat sub-sampling procedure
  
  # perform calculations
  all_ari <- c()
  for (i in 1:nreps) {
    
    # Downsample each individual pc matrix
    ind_downsampled_pcs <- purrr::map(ind_pcs, \(pcs) downsample_pcs_for_metrics(pcs, frac_cells, num_pcs)) 
    
    # downsample integrated pcs so that equivalent fraction of cells from each batch is represented
    # this is required so that the final cluster assignments contain the equivalent number of cells as 
    # the original batch
    integrated_downsampled_pcs <- integrated_pcs |>
      as.data.frame() |>
      tibble::rownames_to_column("cell_id") |> 
      dplyr::mutate(batch = stringr::word(cell_id, -1, sep = "-")) |> 
      dplyr::group_by(batch) |>
      dplyr::slice_sample(prop = frac_cells) |> 
      dplyr::ungroup() |> 
      tibble::column_to_rownames("cell_id") |> 
      dplyr::select(-"batch") |> 
      as.matrix()
    
  
    for (k in k_range) {
      # only cluster integrated pcs once for each k
      integrated_clustering_result <- kmeans(x = integrated_downsampled_pcs, centers=k)
      
      # repeat clustering for each individual pc matrix
      all_ari <- c(all_ari, 
                   purrr::imap_chr(ind_downsampled_pcs, 
                                 \(downsampled_pcs, batch_id){
                                   
                                   # grab clusters for cells from specified batch id
                                   clusters_to_keep <- grep(batch_id, names(integrated_clustering_result$cluster))
                                   batch_integrated_clusters <- unname(integrated_clustering_result$cluster[clusters_to_keep])
                                   
                                   # Perform k-means clustering for individual pcs
                                   ind_clustering_result <- kmeans(x = downsampled_pcs$pcs, centers=k)
                                   ind_clusters <- unname(ind_clustering_result$cluster)
                                   
                                   # Calculate reverse ARI
                                   ari <- pdfCluster::adj.rand.index(batch_integrated_clusters, ind_clusters)
                                   
                                 }))
      
    }
  }
  
  # grab length of batch ids
  batch_id_length <- length(unique(names(all_ari)))
  
  # Create tibble with ari, batch id, and corresponding information
  ari_tibble <- tibble::tibble(
    ari = all_ari,
    batch_id = names(all_ari),
    rep = rep(1:nreps, each=length(k_range)*batch_id_length),
    k   = rep(k_range, nreps*batch_id_length),
    integration_method = integration_method
  )
  
  # Return tibble with batch ARI results which can further be summarized downstream
  return(ari_tibble)
}