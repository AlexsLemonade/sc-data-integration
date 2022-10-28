
# Load libraries function(s) may need -----------
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(magrittr)
})



#' Find all nested files in a given downstream analyses output directory. If 
#'   some files are missing, this function will throw an error and prompt users to
#'   run downstream analyses.
#'
#' @param sce_dir downstream analyses output directory
#' @param library_ids Vector of nested result directory names which corresponds
#'  to library IDs
#'
#' @return Vector of all downstream analyses outputted SCE files, unless an error is thrown.
find_downstream_sce_files <- function(library_ids, sce_dir) {

  # find SCE files that match library ID
  library_search <- paste(library_ids, collapse = "|")
  all_library_files <- list.files(sce_dir,
                                  pattern = library_search,
                                  recursive = TRUE,
                                  full.names = TRUE)
  # just include RDS files, otherwise HTML files will also be found
  sce_files <- all_library_files[grep(pattern = ".rds", all_library_files, ignore.case = TRUE)]
  
  # if the number of sce files is different then the library ID's find the missing files
  if(length(sce_files) < length(library_ids)){
    
    libraries_found <- stringr::str_extract(sce_files, library_search)
    missing_libraries <- setdiff(library_ids, libraries_found)
    
    stop(
      glue::glue(
        "\nMissing SCE object for {missing_libraries}.
      Make sure that you have run `01-run-downstream-analyses.sh` or provided the correct input directory."
      )
    )
  }
  
  
  # Return sce_files
  return(sce_files)
}










#' Combine two or more SCE objects
#'
#' This function combines one or more SCE objects into a single SCE object, with
#'   a cell column indicating the different originating SCEs.
#'
#' @param sce_list Named list of SCE objects to combine, where names are library
#'   biospecimen IDs. No specific assays or dimReduced are expected.
#' @param batch_column Name of column to create that specifies batches. Default
#'   is "batch".
#' @param preserve_rowdata_columns An array of columns that appear in SCE rowData
#'    which should not be renamed with the given library ID. For example,
#'    such an array might be: `c("Gene", "ensembl_ids", "gene_names")`
#'
#' @return The combined SCE object
#'
#' @examples
#' combine_sce_objects(sce_list, c("Gene", "ensembl_ids", "gene_names"))
combine_sce_objects <- function(sce_list = list(),
                                batch_column = "batch",
                                preserve_rowdata_columns = NULL) {


  # Ensure `sce_list` is named (according to library IDs) ----------------------
  if (is.null(names(sce_list))) {
    stop("The `sce_list` must be named by the SCE object's library IDs.")
  }
  # Ensure `sce_list` has >=2 items --------------------------------------------
  if (length(sce_list) < 2) {
    stop("The `sce_list` must contain 2 or more SCE objects to combine.")
  }

  # Ensure that colnames of colData match across all SCE objects ---------------
  # This is required for cbind()
  # Find all the existing column names across SCE objects
  sce_colnames <- sce_list %>%
    purrr::map(~ colnames(colData(.))) %>%
    unname() %>%
    unlist() %>%
    unique()
  # For any SCE object that is *missing* a column,
  #  add in a "dummy" column (NA) of that name
  for (i in 1:length(sce_list)) {
    this_sce_colnames <- colnames(colData(sce_list[[i]]))
    missing_colnames <- setdiff(sce_colnames, this_sce_colnames)
    for (add_col in missing_colnames) {
      colData(sce_list[[i]])[,add_col] <- NA
    }
  }

  # Check that colData colnames now match by comparing all to the first
  new_sce_colnames <- sce_list %>%
    purrr::map(~ colnames(colData(.)))
  for (i in 2:length(sce_list)) {
    # There should be no difference:
    total_diffs <- length(setdiff(new_sce_colnames[[1]], new_sce_colnames[[i]]))
    if (total_diffs != 0) {
      stop("colData column names must match across SCE objects.")
    }
  }

  # Subset SCEs to shared genes and ensure appropriate naming ------------------

  # First, obtain intersection among all SCE objects
  shared_genes <- sce_list %>%
    purrr::map(rownames) %>%
    purrr::reduce(intersect)

  # Now, loop over SCEs to subset each to the array of `shared_genes`
  #  At the same time, we also update the rowData column names to be unique across SCEs
  for (i in 1:length(sce_list)){
    
    # remove any alternative experiments before merging for now
    alt_names <- altExpNames(sce_list[[i]])
    if(!is.null(alt_names)){
      altExp(sce_list[[i]]) <- NULL
    }

    # Find library name
    library_name <- names(sce_list)[i]

    # Subset to shared genes
    sce_list[[i]] <- sce_list[[i]][shared_genes,]

    # Add relevant library ID to rowData column names
    colnames(rowData(sce_list[[i]])) <- paste0(colnames(rowData(sce_list[[i]])), "-", library_name)

    # Add relevant library ID to colData row names
    colnames(sce_list[[i]]) <- paste0(colnames(sce_list[[i]]), "-", library_name)

    # Add relevant library ID to metadata names if metadata is present
    if(length(metadata(sce_list[[i]])) != 0){
      names(metadata(sce_list[[i]])) <- paste0(names(metadata(sce_list[[i]])), "-", library_name) 
    }

    # Add colData column to track this batch
    sce_list[[i]][[batch_column]] <- library_name
  }

  # Combine SCE objects with `cbind()` -----------------------------------------
  combined_sce <- do.call(cbind, sce_list)


  # Restore rowData colnames that do contain shared gene info (not stats) ------
  # In addition, keep only 1 column of these "duplicates"
  for (restore_colname in preserve_rowdata_columns) {

    # Determine which columns need to be updated
    columns <- names(rowData(combined_sce))
    restore_cols <- which(stringr::str_starts(paste0(columns, "-"), restore_colname))

    # Rename the relevant column without the -library_id the first time it appears,
    #  and remove the "duplicated" columns entirely
    names(rowData(combined_sce))[restore_cols[1]] <- restore_colname
    rowData(combined_sce)[restore_cols[-1]] <- NULL
  }


  # Retain only colData names that are the `batch_column` or columns added by
  #  scuttle::addPerCellQC() during the scpca-downstream-analyses processing
  mito_names <- names(colData(combined_sce)) %>%
      stringr::str_subset("^subsets_mito")

  retain_cols <-  c(batch_column,
                    "celltype",
                    # scuttle::addPerCellQC() columns
                    "sum", "detected", mito_names)
  retain_pos <- which(names(colData(combined_sce)) %in% retain_cols)
  colData(combined_sce) <- colData(combined_sce)[, retain_pos]


  # Return combined SCE object ----------------------------
  return(combined_sce)
}


#' Identify highly variable genes for combined SingleCellExperiment objects
#'
#' @param combined_sce SingleCellExperiment object containing normalized gene expression data
#'   from more than one library.
#' @param num_genes Number of highly variable genes to select. Default is 5000.
#' @param batch_column Column present in colData of the SingleCellExperiment object
#'   that contains the original identity of each library. Default is "batch".
#'
#' @return Highly variable gene vector
#'
perform_hvg_selection <- function(combined_sce,
                                  num_genes = 5000,
                                  batch_column = "batch"){

  # check that logcounts are present in combined_sce, required for modelGeneVar
  if(!"logcounts" %in% assayNames(combined_sce)){
    stop("log-normalized counts are not found in the 'logcounts' assay of the combined SCE.")
  }

  # check to make sure that the column to use for blocking is present
  if(!batch_column %in% colnames(colData(combined_sce))){
    stop("batch_column must be a column present in the colData of the SCE object.")
  }

  # extract the column with the block variable
  batch_column <- colData(combined_sce)[,batch_column]

  # model gene variance
  gene_var_block <- scran::modelGeneVar(combined_sce,
                                        block = batch_column)
  # identify subset of variable genes
  gene_list <- scran::getTopHVGs(gene_var_block,
                                 n = num_genes)

  return(gene_list)

}


#' Calculate PCA and UMAP embeddings for a combined SingleCellExperiment object
#'
#' @param combined_sce SingleCellExperiment object containing gene expression data
#'   from more than one library.
#' @param var_genes List of highly variable genes to use for PCA calculation.
#'   Required if PCA `do_pca` is set to TRUE.
#' @param prefix Prefix to use for naming the PCA and UMAP embeddings in the
#'   format `"<prefix>_PCA"` and `"<prefix>_UMAP"` respectively,  that will be
#'   stored in `reducedDim(combined_sce)`. If no prefix is provided, results will
#'   be stored to the `PCA` and `UMAP` slots.
#' @param assay Name of the Assay holding the gene expression matrix to use for
#'   performing PCA. Default is "logcounts". An assay name is only required if
#'   using `pca_type = single`. For `pca_type = multi` the `logcounts` assay is used.
#' @param pca_type Type of PCA to perform prior to UMAP, "single" uses `scater::runPCA()`,
#'   while "multi" uses `batchelor::multiBatchPCA()`. If a PCA method is not chosen,
#'   PCA is skipped, and the existing PCA results must be saved in the object with `prefix`_PCA
#'   and the `prefix` must be used indicated using the `prefix` argument.
#' @param batch_column Column present in colData of the SingleCellExperiment object
#'   that contains the original identity of each library. Default is "batch".
#'
#' @return Combined SingleCellExperiment with PCA and UMAP stored in reducedDim
#'
perform_dim_reduction <- function(combined_sce,
                                  var_genes = NULL,
                                  prefix = NULL,
                                  assay = "logcounts",
                                  pca_type = NULL,
                                  batch_column = "batch"){

  # create pca and umap names
  pca_name <- "PCA"
  umap_name <- "UMAP"
  if(!is.null(prefix)){
    pca_name <- paste(prefix, pca_name, sep = "_")
    umap_name <- paste(prefix, umap_name, sep = "_")
  }


  # check that pca_type is either single or multi if input
  if (!is.null(pca_type)){

    if(!(pca_type %in% c("single", "multi"))){
      stop("pca_type must either be `single` or `multi`.")
    }

    # check that var_genes was provided
    if(is.null(var_genes)){
      stop("A list of variable genes to perform PCA must be provided
           using the var_genes argument.")
    }

    if(pca_type == "multi"){

      # check that logcounts are present in combined_sce, required for modelGeneVar
      if(!"logcounts" %in% assayNames(combined_sce)){
        stop("log-normalized counts are required for multiBatchPCA,
           and are not found in the 'logcounts' assay of the combined SCE.")
      }

      # check to make sure that the column to use for blocking is present
      if(!batch_column %in% colnames(colData(combined_sce))){
        stop("batch_column must be a column present in the colData of the SCE object.")
      }

      # extract the column with the batch_column variable
      batch_column <- colData(combined_sce)[,batch_column]

      # perform multi batch PCA
      multi_pca_list <- batchelor::multiBatchPCA(combined_sce,
                                                 subset.row = var_genes,
                                                 batch = batch_column,
                                                 preserve.single = TRUE)

      # add dataframe with PCA results to SCE object
      reducedDim(combined_sce, pca_name) <- multi_pca_list@listData[[1]]

    }

    if(pca_type == "single"){

      # check for assay present in SCE object
      if(!assay %in% assayNames(combined_sce)){
        stop("assay provided is not in the assayNames() of the SCE object.")
      }

      # add PCA
      combined_sce <- combined_sce %>%
        scater::runPCA(subset_row = var_genes,
                       name = pca_name,
                       exprs_values = assay)

    }

  }

  # add UMAP
  combined_sce <- combined_sce %>%
      scater::runUMAP(dimred = pca_name,
                    name = umap_name)

  return(combined_sce)

}


#' Identify variable genes for a merged object, add to metadata, and return
#' merged object with only gene expression data corresponding to variable genes
#'
#' @param combined_sce SCE object that has been merged using combine_sce_objects
#' @param num_genes Number of highly variable genes to select. Default is set to 5000.
#' @param subset_hvg Indicates whether or not to subset the merged SCE object by highly variable genes.
#'   If --subset_hvg is TRUE, the merged SCE object will only contain genes
#'   identified as highly variable genes.
#' @param batch_column Column present in colData of the SingleCellExperiment object
#'   that contains the original identity of each library. Default is "batch".
#'
#' @return combined SCE object with variable genes added to metadata
set_var_genes <- function(combined_sce,
                          num_genes = 5000,
                          subset_hvg = FALSE,
                          batch_column = "batch"){

  if(!is.logical(subset_hvg)){
    stop("--subset_hvg must be either TRUE or FALSE")
  }

  # grab variable genes
  var_genes <- perform_hvg_selection(combined_sce = combined_sce,
                                     num_genes = num_genes,
                                     batch_column = batch_column)

  # add variable genes to metadata
  metadata(combined_sce)$variable_genes <- var_genes

  if(subset_hvg){
    # subset to only variable genes
    combined_sce <- combined_sce[var_genes,]
  }

  # indicate in metadata if object has been subset by HVG or not
  metadata(combined_sce)$subset_hvg <- subset_hvg

  return(combined_sce)
}


#' Remove original uncorrected expression matrices from an SCE object
#'
#' @param sce_object The SCE object to remove assays from by setting them to NULL
#' @param assays_to_remove Vector of expression assays that will be set to NULL.
#' By default, this includes `counts` and `logcounts`
#'
#' @return The SCE object with specified assays removed
remove_uncorrected_expression <- function(sce_object,
                                          assays_to_remove = c("counts", "logcounts")) {
  for (assay_name in assays_to_remove) {
    assay(sce_object, assay_name) <- NULL
  }
  return(sce_object)
}


#' Add column with cell type annotations to colData of SCE object
#'
#' @param sce_object The SCE object to add cell type annotations
#' @param celltype_info_df Tibble containing celltype specific information for provided SCE.
#'   Columns labeled `barcode` and `celltype` must be present.
#'
#' @return Modified SCE object with cell type annotations added as a column in the colData
#'
add_celltype_info <- function(sce_object,
                              celltype_info_df){
  
  # check that provided celltype info df contains specified columns 
  if(!all(c("barcode", "celltype") %in% colnames(celltype_info_df))){
    stop("`celltype_info_df` must contain the specified columns, `barcode` and `celltype`.")
  }
  
  # check that the barcode column contains overlap with the barcodes in the SCE object
  # right now check that at least some of the barcodes overlap
  if(!any(colnames(sce_object) %in% celltype_info_df$barcode)){
    stop("barcodes column of celltype_info_df provided does not contain barcodes found 
         in the provided SCE object.")
  }
  
  # create new coldata df containing added celltype column
  coldata_df <- as.data.frame(colData(sce_object)) %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(celltype_info_df, by = c("barcode")) %>%
    tibble::column_to_rownames("barcode")
  
  # add coldata back to sce object
  colData(sce_object) <- DataFrame(coldata_df)

  return(sce_object)
}


#' Checks that a given integration method is acceptable
#'
#' @param integration_method The provided integration method to check
#'
#' @return An all-lowercase version of the given provided_method
check_integration_method <- function(integration_method) {
  
  all_integration_methods <- c("fastMNN", "harmony", "rpca", "cca", "scvi", "scanorama")
  if (is.null(integration_method)) {
    stop("An `integration_method` must be provided.")
  } else {
    integration_method <- tolower(integration_method)
    if (!(integration_method %in% tolower(all_integration_methods))) {
      stop(
        paste("The `integration_method` must be one of: ",
              paste(all_integration_methods, collapse = ", ")
        )
      )
    }
  }
  return(integration_method)
}



#' Determines the appropriate name in an SCE object's reduced dimensions slot 
#'  for extracting PCs or analogous reductions returned by other integration methods
#'  
#' @param integration_method The integration method
#'
#' @return The name for the reduced dimension to use
get_reduced_dim_name <- function(integration_method) {
  
  if (integration_method == "scanorama") {
    reduced_dim_name <- "scanorama_SVD"
  } else if (integration_method == "scvi") {
    reduced_dim_name <- "scvi_latent"
  } else {
    reduced_dim_name <- paste0(integration_method, "_PCA")
  }
  
  return(reduced_dim_name)
}



#' Downsample PCs for use in integration metric calculations 
#'  
#' @param integrated_pcs The full set of integrated pcs
#' @param frac_cells The fraction of cells to downsample to
#' @param num_pcs The number of PCs to downsample to
#'
#' @return List with two items: `pcs`, the downsampled PCs; `batch_labels`, the 
#'  corresponding batch labels as integers for the downsampled PCs
downsample_pcs_for_metrics <- function(integrated_pcs, frac_cells, num_pcs) {
  
  num_cells <- nrow(integrated_pcs)
  
  # Determines rows to sample
  downsampled_indices <- sample(1:num_cells,
                                frac_cells*num_cells,
                                replace = FALSE) 
  
  # Extract PCs for downsample, considering only the top `num_pcs`
  downsampled_integrated_pcs <- integrated_pcs[downsampled_indices,1:num_pcs]
  
  
  return (
    list(
      pcs = downsampled_integrated_pcs,
      batch_labels = rownames(downsampled_integrated_pcs)
    )
  )
  
}


#' Remove NAs batches from PCs for use in integration metric calculations 
#'  
#' @param integrated_pcs The full set of integrated pcs
#' @param batches Vector of cell-wise batch information whose length is the number of 
#'   rows in `intergated_pcs`
#' 
#' @return List with two items: `pcs`: PCs with NA batch cells removed and 
#'   rownames assigned as batch; `indices`: Which rows were _retained_ 
remove_batch_nas_from_pcs <- function(integrated_pcs, batches) {
  
  retain_indices <- which(!is.na(batches))
  batches <- batches[retain_indices]
  integrated_pcs <- integrated_pcs[retain_indices,]
  
  # check dimensions still match:
  if (nrow(integrated_pcs) != length(batches)) {
    stop("Incompatable PC and batch information dimensions after removing NAs.")
  }
  
  # Set PC rownames to be the batches, for use in obtaining downsampled labels
  rownames(integrated_pcs) <- batches
  
  return(
    list(
      pcs = integrated_pcs,
      indices = retain_indices
    )
  )
}
  
