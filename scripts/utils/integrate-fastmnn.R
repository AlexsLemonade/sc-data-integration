library(SingleCellExperiment) # Needed for assays() function


#' Integrated combined SCE objects with `fastMNN`
#'
#' This function integrates a combined SCE object the using `fastMNN` function
#'  from the `batchelor` package.
#' @param combined_sce The combined SCE objects to integrate
#' @param gene_list Vector of high-variance genes to consider. The default value 
#'   of `NULL` means all genes will be used.
#' @param cosine_norm Logical indicating whether cosine normalization should be 
#'   performed prior to calculating PCs for integration. Default: `TRUE`. 
#' @param fastmnn_k Number of nearest-neighbors to consider when identifying 
#' mutual nearest neighbors. Default: 20 (same default as in `batchelor::fastMNN()`)
#' @param fastmnn_d Number of PCs to use when performing PCA. 
#'   Default: 50 (same default as in `batchelor::fastMNN()`)
#' @param seed Random seed to set for `fastMNN` integration
#' @param ... Additional arguments to pass into `batchelor::fastMNN()`
#'
#' @return The integrated SCE object
integrate_fastMNN <- function(combined_sce, 
                              gene_list = NULL, 
                              cosine_norm = TRUE,
                              fastmnn_k = 20, 
                              fastmnn_d = 50, 
                              seed = 2022,
                              ...) {
  
  set.seed(seed)
  
  # Throw an error if normalization has not been performed
  if (!("logcounts" %in% names(assays(combined_sce)))) {
    stop("The `combined_sce` object is missing a `logcounts` assay required for fastMNN integration.")
  }
  
  # Perform integration with fastMNN -------------------
  integrated_sce <- batchelor::fastMNN(combined_sce, 
                                       # Specify batches. TODO: More batches?
                                       batch = combined_sce$batch,
                                       # Which genes to use for integration
                                       # The default value of NULL uses all genes
                                       subset.row = gene_list,
                                       # Perform cosine normalization?
                                       cos.norm = cosine_norm,
                                       # How many nearest neighbors?
                                       k = fastmnn_k,
                                       # How many PCs?
                                       d = fastmnn_d,
                                       # Anything else?
                                       ...)
  
  # Add sce_integrated fields into `combined_sce` ---------------
  
  # From docs: the `reconstructed` assay is a...
  #  low-rank reconstruction of the expression matrix. 
  #  This can be interpreted as per-gene corrected log-expression values 
  #  (after cosine normalization, if cos.norm=TRUE) but should not be 
  #  used for quantitative analyses.
  assay(combined_sce, "fastMNN_corrected")  <- assay(integrated_sce, "reconstructed")
  # The integrated PCs
  reducedDim(combined_sce, "fastMNN_PCA") <- reducedDim(integrated_sce, "corrected")
  
  # Return SCE object with fastMNN information 
  return(combined_sce)
  
}
