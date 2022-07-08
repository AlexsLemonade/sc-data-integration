library(SingleCellExperiment) # Needed for assays() function



#' integrate_fastMNN
#'
#' Integrate SCE datasets using fastMNN from the batchelor package
#' @param sce1 The first SCE object to integrate
#' @param sce2 The second SCE object to integrate
#' @param gene_list Optional vector of high-variance genes to integrate on. Default: NULL (all genes will be used)
#' @param use_all_genes Optional logical indicating whether all genes should be used. Default: FALSE. This argument is 
#' @param num_genes Optional number of high-variance genes to identify for use in integration. Default: 5000. This argument is *ignored* if gene_list is not NULL and/or if use_all_genes is TRUE 
#' @param fastmnn_k Number of nearest-neighbors to consider when identifying mutual nearest neighbors. Default: 20 (same default as in batchelor::fastMNN())
#' @param fastmnn_d Number of PCs to use when performing PCA. Default: 50 (same default as in batchelor::fastMNN())
#' @param seed Random seed to set for integration
#' @param ... Additional arguments to pass into batchelor::fastMNN()
#'
#' @return List containing the integrated SCE object (`sce_integrated`) and the unintegrated SCE objects (`sce1` and `sce2`)
#'
#' @examples
integrate_fastMNN <- function(sce1, sce2, 
                              gene_list = NULL, 
                              num_genes = 5000, 
                              use_all_genes = FALSE,
                              fastmnn_k = 20, 
                              fastmnn_d = 50, 
                              seed = 2022,
                              ...) {
  
  set.seed(seed)
  
  # Add logcounts assay to SCE objects if it does not exist
  if (!("logcounts" %in% names(assays(sce1)))) {
    sce1 <- scuttle::logNormCounts(sce1)
  }
  if (!("logcounts" %in% names(assays(sce2)))) {
    sce2 <- scuttle::logNormCounts(sce2)
  }
  common_genes <- intersect(rownames(sce1), rownames(sce2))
  common_sce1 <- sce1[common_genes,]
  common_sce2 <- sce2[common_genes,]
  # Set up gene list:
  #   If gene_list is not NULL, the user-given vector will be used (no checking is performed here!)
  #   If gene_list is NULL....
  #      And use_all_genes is TRUE, all genes are used
  #      And use_all_genes is FALSE, `num_genes` high variance genes are used (via scran)
  
  # If no genes provided *and* we don't want to use all genes, 
  #  allow scran to identify high-variance genes to use based on provided `num_genes`
  if (is.null(gene_list) & !(use_all_genes)) {
    # Use vignette procedure to pull out high-variance genes 
    #  with user-given num_genes, or a default of 5000
    dec1 <- scran::modelGeneVar(sce1)
    dec2 <- scran::modelGeneVar(sce2)
    combined.dec <- scran::combineVar(dec1, dec2)
    gene_list <- scran::getTopHVGs(combined.dec, n = num_genes)   
  } 
  
  # Perform integration with fastMNN
  integrated <- batchelor::fastMNN(sce1, sce2, 
                                   # Which genes to use for integration (NULL uses all genes)
                                   subset.row = gene_list,
                                   # How many nearest neighbors?
                                   k = fastmnn_k,
                                   # How many PCs?
                                   d = fastmnn_d,
                                   # Anything else?
                                   ...)
  integrated 
  
  # Return integrated SCE object as well as 
  # the non-integrated SCE objects that now have logcounts
  return( 
    list(
      sce_integrated = integrated,
      sce1 = sce1,
      sce2 = sce2
    ))
  
}
