#' Integrated combined SCE objects with Seurat's `IntegrateData`
#' 
library(SingleCellExperiment)
library(Seurat)
library(magrittr)

#' Integrate seurat objects with Seurat reciprocal PCA or CCA
#'
#' @param seurat_list List of seurat objects normalized with Seurat::SCTransform()
#' @param reduction_method The Seurat reduction method to implement with 
#' Seurat::FindIntegrationAnchors(), can be "rpca" or "cca"
#' @param num_genes Number of variable features to use for integration
#' @param batch_column The column variable indicating batches, which would 
#' typically corresponds to the library ID. Default is "batch".
#' @param integration_dims The range of values to supply to the `dims` argument
#' when implementing the integration functions
#' @param umap_dims The range of values to supply to the `dims` argument when
#' implementing the `runUMAP()` function. Default is 1:10.
#'
#' @return integrated seurat object
#'
#'

integrate_seurat <- function(seurat_list,
                             reduction_method,
                             num_genes = 3000,
                             batch_column = "batch",
                             integration_dims,
                             umap_dims = c(1:10),
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
  anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                            reduction = reduction_method,
                                            dims = integration_dims)
  
  
  integrated_object <- Seurat::IntegrateData(anchorset = anchors,
                                             normalization.method = "SCT",
                                             features.to.integrate = shared_genes,
                                             dims = integration_dims)
  
  # set idents as the batch column
  Seurat::Idents(integrated_object) <- batch_column
  
  # make sure that the default assay is set to integrated
  Seurat::DefaultAssay(integrated_object) <- "integrated"
  
  # add PCA and UMAP to the integrated object
  integrated_object <- integrated_object %>%
    Seurat::RunPCA(reduction.name = paste0("seurat", reduction_method, "_PCA")) %>%
    Seurat::RunUMAP(reduction.name = paste0("seurat", reduction_method, "_UMAP"),
                    reduction = paste0("seurat", reduction_method, "_PCA"), dims = umap_dims)
  
  return(integrated_object)
  
}
