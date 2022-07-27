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
#' @param num_features Number of variable features to use for integration
#'
#' @return integrated seurat object
#'
#'

integrate_seurat <- function(seurat_list,
                             reduction_method,
                             num_features = 3000){
  
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
  
  # find common variable features for integration
  common_features <- Seurat::SelectIntegrationFeatures(seurat_list, nfeatures = num_features)
  
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
                                            dims = 1:50)
  
  
  integrated_object <- Seurat::IntegrateData(anchorset = anchors,
                                             normalization.method = "SCT",
                                             features.to.integrate = shared_genes,
                                             dims = 1:50)
  
  # set idents as library ID
  Seurat::Idents(integrated_object) <- "library"
  
  # make sure that the default assay is set to integrated
  Seurat::DefaultAssay(integrated_object) <- "integrated"
  
  # add PCA and UMAP to the integrated object
  integrated_object <- integrated_object %>%
    Seurat::RunPCA() %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:10)
  
  return(integrated_object)
  
}
