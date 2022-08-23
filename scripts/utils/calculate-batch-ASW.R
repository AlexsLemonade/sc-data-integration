suppressPackageStartupMessages({
  library(magrittr)
  library(bluster)
  library(SingleCellExperiment)
})

#' Function to calculate batch ASW scores from an integrated SCE object
#'
#' This function uses a similar approach as used in this paper:
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
#' Data is subsampled to 80% of cells. A given number of PCs (default 20) are
#' used for distance calculations to obtain ASW scores.
#' This procedure is repeated 20x.
#'
#' @param integrated_sce The integrated SCE object
#' @param num_pcs The number of PCs to use during k-means clustering. Default: 20
#' @param seed Seed for initializing random sampling
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#'
#' @return Tibble with three columns: `rep`, representing the given downsampling replicate;
#'   `k`, the given k for k-means; `batch_ari`, the calculated batch ARI for the given
#'   combination of `rep` and `k`
calculate_batch_asw <- function(integrated_sce,
                                num_pcs = 20,
                                seed = NULL,
                                integration_method = NULL) {

  # Check integration method
  integration_method <- check_integration_method(integration_method)

  # Set the seed for subsampling 
  set.seed(seed)

  # Set up parameters
  frac_cells <- 0.8        # fraction of cells to downsample to
  nreps <- 20              # number of times to repeat sub-sampling procedure

  # pull out relevant information for calculations
  reduced_dim_name <- get_reduced_dim_name(integration_method)
  all_integrated_pcs <- reducedDim(integrated_sce, reduced_dim_name)

  # perform calculations
  all_batch_asw <- tibble::tibble(
    rep         = as.numeric(),
    cell_name   = as.character(),
    batch_asw   = as.numeric(),
    asw_cluster = as.numeric(),
    asw_other   = as.numeric()
  )
  for (i in 1:nreps) {

    # Downsample PCs
    downsampled <- downsample_pcs_for_metrics(all_integrated_pcs, frac_cells, num_pcs)

    # Calculate batch ASW and add into final tibble
    batch_asw <- bluster::approxSilhouette(downsampled$pcs, downsampled$batch_labels) %>%
      tibble::as_tibble(rownames = "cell_name") %>%
      dplyr::mutate(rep = i) %>%
      dplyr::select(rep, 
                    cell_name, 
                    batch_asw = width, 
                    asw_cluster = cluster, 
                    asw_other = other)
    
    all_batch_asw <- dplyr::bind_rows(all_batch_asw, batch_asw)
  }
    

  # Return tibble with batch ASW results which can further be summarized downstream
  return(all_batch_asw)

}
