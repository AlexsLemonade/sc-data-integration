suppressPackageStartupMessages({
  library(magrittr)
  library(bluster)
  library(SingleCellExperiment)
})

#' Function to calculate silhouette width scores from an integrated SCE object
#'
#' This function uses a similar approach as used in this paper:
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
#' Data is subsampled to 80% of cells. A given number of PCs (default 20) are
#' used for distance calculations to obtain silhouette widths.
#' This procedure is repeated 20x.
#'
#' @param integrated_sce The integrated SCE object
#' @param num_pcs The number of PCs to use during k-means clustering. Default: 20
#' @param seed Seed for initializing random sampling
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param unintegrated Indicates whether the provided data is integrated (`FALSE`; default) or
#'   integrated (`TRUE`).
#' @param batch_column The variable in `integrated_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "batch".
#'   
#' @return Tibble with six columns: `rep`, representing the given downsampling replicate;
#'   `silhouette_width`, the calculated silhouette width for the given `rep`; `silhouette_cluster`, 
#'   the assigned cluster for the cell during silhouette calculation, i.e. the true identity; 
#'   `other_cluster`, the other assigned for the cell during silhouette calculation; 
#'   `integration_method`, the given integration method
calculate_silhouette_width <- function(integrated_sce,
                                       num_pcs = 20,
                                       seed = NULL,
                                       integration_method = NULL, 
                                       unintegrated = FALSE, 
                                       batch_column = "batch") {

  # Set the seed for subsampling 
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
  all_silhouette <- tibble::tibble(
    rep                = as.numeric(),
    silhouette_width   = as.numeric(),
    silhouette_cluster = as.character(),
    other_cluster      = as.character()
  )
  for (i in 1:nreps) {
    # Downsample PCs
    downsampled <- downsample_pcs_for_metrics(final_pcs, frac_cells, num_pcs)

    # Calculate batch ASW and add into final tibble
    rep_silhouette <- bluster::approxSilhouette(downsampled$pcs, downsampled$batch_labels) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(rep = i) %>%
      dplyr::select(rep, 
                    silhouette_width = width, 
                    silhouette_cluster = cluster, 
                    other_cluster = other)
    
    all_silhouette <- dplyr::bind_rows(all_silhouette, rep_silhouette)
  }
    
  # Add integration method into tibble
  all_silhouette <- dplyr::mutate(all_silhouette, integration_method = integration_method)

  # Return tibble with silhouette width results which can further be summarized downstream
  return(all_silhouette)

}
