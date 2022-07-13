library(SingleCellExperiment) # Needed for assays() function


combined_sce <- readr::read_rds(here::here("combined_sce.RDS"))

#' integrate_fastMNN
#'
#' Integrate SCE datasets using fastMNN from the batchelor package
#' @param combined_sce The combined SCE objects to integrate
#' @param use_all_genes Logical indicating whether all genes should be used. Default: FALSE.
#' @param num_genes Optional number of high-variance genes to identify for use in integration. Default: 5000. 
#'   This argument is *ignored* if gene_list if use_all_genes is TRUE 
#' @param fastmnn_k Number of nearest-neighbors to consider when identifying mutual nearest neighbors. Default: 20 (same default as in batchelor::fastMNN())
#' @param fastmnn_d Number of PCs to use when performing PCA. Default: 50 (same default as in batchelor::fastMNN())
#' @param seed Random seed to set for fastMNN integration
#' @param ... Additional arguments to pass into `batchelor::fastMNN()`
#'
#' @return The integrated SCE object
#'
#' @examples
integrate_fastMNN <- function(combined_sce, 
                              num_genes = 5000, 
                              use_all_genes = FALSE,
                              fastmnn_k = 20, 
                              fastmnn_d = 50, 
                              seed = 2022,
                              ...) {
  
  set.seed(seed)
  
  # Throw an error if normalization has not been performed
  if (!("logcounts" %in% names(assays(combined_sce)))) {
    stop("The `combined_sce` object is missing a `logcounts` assay which is required for fastMNN integration.")
  }

  # Set up gene list:
  if (!(use_all_genes)) {
    # Determine high-variance genes if we are not using them all
    gene_var <- scran::modelGeneVar(combined_sce)
    gene_list <- scran::getTopHVGs(gene_var, n = num_genes)   
  } else {
    # Set gene_list to NULL to use all genes, consistent with fastMNN function
    gene_list <- NULL 
  }
  
  # Perform integration with fastMNN
  sce_integrated <- batchelor::fastMNN(combined_sce, 
                                   # Which genes to use for integration (NULL uses all genes)
                                   subset.row = gene_list,
                                   # How many nearest neighbors?
                                   k = fastmnn_k,
                                   # How many PCs?
                                   d = fastmnn_d,
                                   # Anything else?
                                   ...)
  
  # Return integrated SCE object
  return(sce_integrated)
  
}
