library(SingleCellExperiment) 

source(
  here::here(
    "scripts", 
    "utils",
    "integration-helpers.R"
  )
)


#' Integrated combined SCE objects with `fastMNN`
#'
#' This function integrates a combined SCE object the using `fastMNN` function
#'  from the `batchelor` package.
#' @param combined_sce The combined SCE objects to integrate
#' @param batch_column The variable in `combined_sce` indicating batches, typically corresponds to the library ID. Default 
#'   is "batch".
#' @param gene_list Vector of high-variance genes to consider. The default value 
#'   of `NULL` means all genes will be used.
#' @param cosine_norm Boolean indicating whether cosine normalization should be 
#'   performed prior to calculating PCs for integration. Default: `TRUE`. 
#' @param seed Random seed to set for `fastMNN` integration. A seed will only
#'  be set if this is not `NULL` (the default).
#' @param return_uncorrected_expression Boolean indicating whether to return the
#'  uncorrected expression values in the returned SCE object. Default: `FALSE`.
#'  TODO: This is the opposite of what Ally implemented for scanorama, so make sure the final logic is consistent!!
#' @param ... Additional arguments to pass into `batchelor::fastMNN()`
#'
#' @return The integrated SCE object
integrate_fastMNN <- function(combined_sce, 
                              batch_column = "batch",
                              gene_list = NULL, 
                              cosine_norm = TRUE,
                              seed = NULL,
                              return_uncorrected_expression = FALSE,
                              ...) {
  
  if (!(is.null(seed))) {
    set.seed(seed)
  }
  
  # Perform checks ---------------
  if (!("logcounts" %in% names(assays(combined_sce)))) {
    stop("The `combined_sce` object requires a `logcounts` assay for fastMNN integration.")
  }
  
  if (!(batch_column %in% names(colData(combined_sce)))) {
    stop("The provided `batch_column` column must be in `combined_sce` colData.")
  }


  # Perform integration with fastMNN -------------------
  integrated_sce <- batchelor::fastMNN(combined_sce, 
                                       # Specify batches.
                                       batch = colData(combined_sce)[,batch_column],
                                       # Which genes to use for integration
                                       # The default value of NULL uses all genes
                                       subset.row = gene_list,
                                       # Perform cosine normalization?
                                       cos.norm = cosine_norm,
                                       # Anything else?
                                       ...)
  
  # Add integrated_sce fields into `combined_sce` ---------------
  
  # From fastMNN docs: the `reconstructed` assay is a...
  #  low-rank reconstruction of the expression matrix. 
  #  This can be interpreted as per-gene corrected log-expression values 
  #  (after cosine normalization, if cos.norm=TRUE) but should not be 
  #  used for quantitative analyses.

  # We will use `_reconstructed` for fastMNN's `reconstructed` and 
  #  `_PCA` for fastMNN's `corrected`.
  
  # Only add this information if all genes were used, meaning dimensions are compatible
  if (is.null(gene_list)){
    assay(combined_sce, "fastMNN_reconstructed")  <- assay(integrated_sce, "reconstructed")
  }
  
  # Add in the PCs, regardless of the gene list
  reducedDim(combined_sce, "fastMNN_PCA") <- reducedDim(integrated_sce, "corrected")
  
  # Perform UMAP with the new PCs -----------------
  combined_sce <- perform_dim_reduction(combined_sce, prefix = "fastMNN")
  
  # Remove uncorrected expression values, unless otherwise specified ----
  if (!return_uncorrected_expression) {
    assay(combined_sce, "counts") <- NULL
    assay(combined_sce, "logcounts") <- NULL
  }
  
  # Return SCE object with fastMNN information ---------------
  return(combined_sce)
  
}
