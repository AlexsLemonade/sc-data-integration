
# Load libraries function(s) may need -----------
library(SingleCellExperiment)
library(magrittr)


#' Combine two or more SCE objects
#' 
#' This function combines one or more SCE objects into a single SCE object, with
#'   a cell column `batch` indicating the different originating SCEs.
#'
#' @param sce_list Named list of SCE objects to combine, where names are library 
#'   biospecimen IDs. No specific assays or dimReduced are expected.
#' @param preserve_rowdata_columns An array of columns that appear in SCE rowData
#'    which should not be renamed with the given library ID. For example, 
#'    such an array might be: `c("Gene", "ensembl_ids", "gene_names")`
#'
#' @return The combined SCE object
#'
#' @examples 
#' combine_sce_objects(sce_list, c("Gene", "ensembl_ids", "gene_names"))
combine_sce_objects <- function(sce_list = list(), 
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
    
    # Subset to shared genes
    sce_list[[i]] <- sce_list[[i]][shared_genes,]
      
    # Add relevant sample IDs to rowData column names
    colnames(rowData(sce_list[[i]])) <- paste0(colnames(rowData(sce_list[[i]])), "-", library_ids[i])
    
    # Add relevant sample IDs to colData row names
    colnames(sce_list[[i]]) <- paste0(colnames(sce_list[[i]]), "-", library_ids[i])
    
    # Add colData column to track this batch
    sce_list[[i]]$batch <- library_ids[i]
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

  

  # Return combined SCE object ----------------------------
  return(combined_sce)
}


#' Identify highly variable genes for combined SingleCellExperiment objects
#'
#' @param combined_sce SingleCellExperiment object containing normalized gene expression data 
#'   from more than one library.
#' @param num_genes Number of highly variable genes to select. Default is 5000.
#' @param block_var Column present in colData of the SingleCellExperiment object 
#'   that contains the original identity of each library. Default is "batch". 
#'
#' @return Highly variable gene vector
#'
perform_hvg_selection <- function(combined_sce,
                                  num_genes = 5000,
                                  block_var = "batch"){
  
  # check that logcounts are present in combined_sce, required for modelGeneVar
  if(!"logcounts" %in% assayNames(combined_sce)){
    stop("log-normalized counts are not found in the 'logcounts' assay of the combined SCE.")
  }
  
  # check to make sure that the column to use for blocking is present 
  if(!block_var %in% colnames(colData(combined_sce))){
    stop("block_var must be a column present in the colData of the SCE object.")
  }
  
  # extract the column with the block variable
  block_col <- colData(combined_sce)[,block_var]
  
  # model gene variance 
  gene_var_block <- scran::modelGeneVar(combined_sce, 
                                        block = block_col)
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
#'   `single_pca` is TRUE. For `multi_pca` the `logcounts` assay is used.
#' @param multi_pca Boolean indicating whether or not to perform multiBatchPCA prior to UMAP.
#'   Default is set to TRUE. This is the recommended method for combined SCE objects 
#'   prior to integration. Note that setting both multi_pca and single_pca to FALSE skips PCA, 
#'   and the existing PCA results must be saved in the object with `prefix`_PCA 
#'   and the `prefix` must be used indicated using the `prefix` argument.
#' @param single_pca Boolean indicating whether or not to perform PCA using `scater::runPCA`
#'   on the combined SCE prior to UMAP. Default is set to FALSE. This is the 
#'   recommneded method for combined SCE objects after integration is complete 
#'   if PCA embeddings have not already been calculated.Note that setting both 
#'   multi_pca and single_pca to FALSE skips PCA, and the existing PCA results 
#'   must be saved in the object with `prefix`_PCA and the `prefix` must be used 
#'   indicated using the `prefix` argument.
#' @param batch_var Column present in colData of the SingleCellExperiment object 
#'   that contains the original identity of each library. Default is "batch". 
#'
#' @return Combined SingleCellExperiment with PCA and UMAP stored in reducedDim
#'
perform_dim_reduction <- function(combined_sce, 
                                  var_genes = NULL, 
                                  prefix = NULL, 
                                  assay = "logcounts",
                                  multi_pca = TRUE,
                                  single_pca = FALSE,
                                  batch_var = "batch"){
  
  # check that only one PCA type is selected 
  if (multi_pca && single_pca){
    stop("Please only specify one type of PCA to perform.
         If using single PCA be sure to turn off multi PCA with multi_pca = FALSE.")
  }
  
  # create pca and umap names 
  pca_name <- "PCA"
  umap_name <- "UMAP"
  if(!is.null(prefix)){
    pca_name <- paste(prefix, pca_name, sep = "_")
    umap_name <- paste(prefix, umap_name, sep = "_")
  }
  
  # if performing either PCA check that variable genes are provided
  if(multi_pca || single_pca){
    # check that var_genes was provided
    if(is.null(var_genes)){
      stop("A list of variable genes to perform PCA must be provided 
           using the var_genes argument.")
    }
  }
  
  if(multi_pca){
    
    # check that logcounts are present in combined_sce, required for modelGeneVar
    if(!"logcounts" %in% assayNames(combined_sce)){
      stop("log-normalized counts are required for multiBatchPCA, 
           and are not found in the 'logcounts' assay of the combined SCE.")
    }
    
    # check to make sure that the column to use for blocking is present 
    if(!batch_var %in% colnames(colData(combined_sce))){
      stop("batch_var must be a column present in the colData of the SCE object.")
    }
    
    # extract the column with the block variable
    batch_col <- colData(combined_sce)[,batch_var]
    
    # perform multi batch PCA 
    multi_pca_list <- batchelor::multiBatchPCA(combined_sce,
                                               subset.row = var_genes,
                                               batch = batch_col,
                                               preserve.single = TRUE)
    
    # add dataframe with PCA results to SCE object
    reducedDim(combined_sce, pca_name) <- multi_pca_list@listData[[1]]
    
  }
  
  if(single_pca){
    
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
  
  # add UMAP 
  combined_sce <- combined_sce %>%
      scater::runUMAP(dimred = pca_name,
                    name = umap_name)
  
  return(combined_sce)
  
}






