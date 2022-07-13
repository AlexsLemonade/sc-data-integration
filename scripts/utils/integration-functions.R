
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







