# Integrated list of Seurat objects with Seurat's `IntegrateData`
# 

suppressPackageStartupMessages({
    library(Seurat)
    library(magrittr)
})



#' Convert a merged SCE object into a list of Seurat objects that are normalized
#'  with `Seurat::SCTransform()`
#'
#' @param combined_sce A merged SCE object containing two or more libraries to 
#'  be integrated
#' @param batch_column The column variable indicating batches, which would 
#' typically corresponds to the library ID. Default is "batch".
#'
#' @return List of normalized Seurat objects to be integrated
#' 
prepare_seurat_list <- function(combined_sce, 
                                batch_column = "batch") {
  
  # Convert and normalize
  seurat_list <- combined_sce %>%
    # convert to Seurat
    scpcaTools::sce_to_seurat() %>%
    # split into list
    # note that the batch_column metadata column is retained during the split
    Seurat::SplitObject(split.by = batch_column) %>%
    # normalize each Seurat object in the list
    purrr::map(Seurat::SCTransform)

  return(seurat_list)
}



#' Integrate seurat objects with Seurat reciprocal PCA or CCA
#'
#' @param seurat_list List of seurat objects normalized with Seurat::SCTransform()
#' @param reduction_method The Seurat reduction method to implement with 
#' Seurat::FindIntegrationAnchors(), can be "rpca" or "cca"
#' @param num_genes Number of variable features to use for integration. Default
#' is 2000, same as Seurat's default.
#' @param batch_column The column variable indicating batches, which would 
#' typically corresponds to the library ID. Default is "batch".
#' @param integration_dims The range of values to supply to the `dims` argument
#' when implementing the integration functions.  Default is 1:30, same as 
#' Seurat's default.
#' @param umap_dims The range of values to supply to the `dims` argument when
#' implementing the runUMAP() function. Default is 1:30, same as implemented
#' in the single cell integration vignette found at 
#' https://www.singlecellcourse.org/scrna-seq-dataset-integration.html#seurat-v3-3-vs-5-10k-pbmc
#' @param anchor_threshold Minimum threshold for number of neighbors to consider when weighting 
#' anchors during integration
#' @param ... Allows for any additional parameters that a user may want to pass
#' to the Seurat::IntegrateData() function
#'
#' @return integrated seurat object with `rpca_` and `cca_` prefixed reducedDims slots
#'
#'
integrate_seurat <- function(seurat_list,
                             reduction_method,
                             num_genes = 2000,
                             batch_column = "batch",
                             integration_dims = 1:30,
                             umap_dims = 1:30,
                             anchor_threshold = 100,
                             ...){
  
  # check that all objects contain the `SCT` assay as default assay
  for(seurat_obj in seurat_list){
    if(Seurat::DefaultAssay(seurat_obj) != "SCT"){
      if("SCT" %in% names(seurat_obj@assays)){
        stop("Normalization must be normalized using Seurat::SCTransform")
      } else {
        # if default assay is not SCT and it's present make sure default is set
        Seurat::DefaultAssay(seurat_obj) <- "SCT"
      }
    }
  }
  
  # check that the batch column exists
  if (!(batch_column %in% names(seurat_obj@meta.data))) {
    stop("The provided `batch_column` column must be in the seurat object's metadata.")
  }
  
  # check that reduction method is valid
  reduction_method <- tolower(reduction_method)
  if (!reduction_method %in% c("cca", "rpca")) {
    stop("The `reduction_method` must be one of `cca` or `rpca` (case-insentitive).")
  }
  
  # find common variable features for integration
  common_features <- Seurat::SelectIntegrationFeatures(seurat_list, nfeatures = num_genes)
  
  # find the shared genes across all objects
  shared_genes <- seurat_list %>%
    purrr::map(rownames) %>%
    purrr::reduce(intersect)
  
  # prepare the normalized data in each object for anchor identification
  seurat_list <- Seurat::PrepSCTIntegration(seurat_list, anchor.features = common_features)
  seurat_list <- seurat_list %>%
    purrr::map(Seurat::RunPCA, features = common_features, verbose = FALSE)
  
  # find anchors and perform integration
  # note that Seurat implements a default seed value of 42 in the applicable
  # functions below, to ensure that the results are reproducible
  anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                            reduction = reduction_method,
                                            dims = integration_dims)
  
  # create a table of anchors found across datasets 
  anchor_table <- table(anchors@anchors[,c("dataset1", "dataset2")])
  # get the minimum number of anchors, taking the second lowest as the lowest will always be 0 across the diagonal
  min_anchors <- unique(sort(anchor_table))[2]
  
  if(min_anchors < anchor_threshold){
    k_weight <- min_anchors
  } else {
    k_weight <- anchor_threshold
  }
  
  integrated_object <- Seurat::IntegrateData(anchorset = anchors,
                                             normalization.method = "SCT",
                                             features.to.integrate = shared_genes,
                                             dims = integration_dims,
                                             k.weight = k_weight,
                                             ...)
  
  # set idents as the batch column
  Seurat::Idents(integrated_object) <- batch_column
  
  # make sure that the default assay is set to integrated
  Seurat::DefaultAssay(integrated_object) <- "integrated"
  
  # add PCA and UMAP to the integrated object
  integrated_object <- integrated_object %>%
    Seurat::RunPCA(reduction.name = paste0(reduction_method, "_PCA")) %>%
    Seurat::RunUMAP(reduction.name = paste0(reduction_method, "_UMAP"),
                    reduction = paste0(reduction_method, "_PCA"), 
                    dims = umap_dims)
  
  return(integrated_object)
  
}
