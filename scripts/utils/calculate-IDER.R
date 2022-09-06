suppressPackageStartupMessages({
  library(CIDER)
  library(SingleCellExperiment)
  library(Seurat)
})

source(
  here::here(
    "scripts",
    "utils",
    "integration-helpers.R"
  )
)

#' Function to calculate the IDER similarity matrix from an integrated SCE object
#'
#' @param integrated_sce The integrated SCE object
#' @param batch_column The variable in `integrated_sce` indicating batches. Default
#'   is "batch".
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param unintegrated Indicates whether the provided data is integrated (`FALSE`; default) or
#'   integrated (`TRUE`).
#'
#' @return Tibble with five columns: `k0_fraction`, the fraction of total sample sized used to
#'  assign the neighborhood size used by `kBET`; `rep`, the given subsample for that `k0`;
#'  `integration_method`,the given integration method; `kbet_stat_type`, either 
#'  "observed_rejection_rate" or "expected_rejection rate"; `kbet_stat`, the given rejection
#'  rate value seen in `kbet_stat_type` for the given subsample/k0 combination
#'  
calculate_IDER <- function(integrated_sce,
                           batch_column = "batch",
                           integration_method = NULL,
                           unintegrated = FALSE,
                           seed = NULL) {

  # Set seed if given
  set.seed(seed)

  # Check integration method, and convert SCE to Seurat
  # Settings depending on whether data is integrated or not
  if (unintegrated){
    # In the end, we'll return "unintegrated" in the data frame for integration method
    integration_method <- "unintegrated"
    
    # use simply "PCA" for reduced dimensions
    reduced_dim_name <- "PCA"
    
    # Use the "logcounts" assay to use in seurat conversion
    assay_name <- "logcounts"
  } else {
    # Check integration method
    integration_method <- check_integration_method(integration_method)
    
    # Get name for reduced dimensions
    reduced_dim_name <- get_reduced_dim_name(integration_method)
    
    # get assay name to use in seurat conversion
    assay_name <- paste0(integration_method, "_corrected")
  }
  
  
  # Convert integrated SCE to Seurat for CIDER input
  seurat_obj <- scpcaTools::sce_to_seurat(integrated_sce, 
                                          assay_name = assay_name)
  # To do: have to convert this mess??
  #> seurat_obj@meta.data
  #[1] orig.ident                 nCount_harmony_corrected   nFeature_harmony_corrected
  #[4] sum                        detected                   subsets_mito_sum          
  #[7] subsets_mito_detected      subsets_mito_percent       celltype                  
  #[10] batch 
  
  # Add PCs into seurat object
  seurat_pcs <- CreateDimReducObject(
    embeddings = reducedDim(integrated_sce, reduced_dim_name),
    assay = assay_name, # must match
    key = "PCA_" # seurat is going to add an underscore here, so let's explicitly include one
  ) 
  
  
  seurat_obj[["pca"]] <- seurat_pcs # creates: seurat_obj@reductions$pca@cell.embeddings
  
  # Ensure "Batch" (case sensitive) is present
  seurat_obj[["Batch"]] <- colData(integrated_sce)[,batch_column]
  
  # Using this approach: 
  # https://github.com/zhiyhu/CIDER#cider-as-an-evaluation-metric---quick-start
  seurat_obj <- hdbscan.seurat(seurat_obj)
  ider <- getIDEr(seu.integrated, verbose = FALSE)
  seu.integrated <- estimateProb(seu.integrated, ider)
  
    
}
