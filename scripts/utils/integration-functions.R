
# Load libraries function(s) may need -----------
library(SingleCellExperiment)
library(magrittr)


#' combine_sce_objects
#'
#' @param sce_list Named list of SCE objects to combine, where names are biospecimen IDs. No specific assays or dimReduced are expected.
#'
#' @return List with two items: `sce_list_updated`, the updated input list with new rowData column names and new colData row names, subsetted to shared genes as specified, and ii) `combined_sce` the SCE objects combined with `cbind()`
#'
#' @examples
combine_sce_objects <- function(sce_list = list()) {
  
  
  # Ensure `sce_list` is named (according to sample IDs) -----------------------
  if (is.null(names(sce_list))) {
    stop("The `sce_list` must be named by the SCE object's sample IDs.")
  }
  
  # Ensure that colnames of colData match across all SCE objects ---------------
  #  which is required for cbind() 
  # Find all the existing column names across SCE objects
  sce_colnames <- sapply(sce_list, function(x) colnames(colData(x))) %>%
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
  new_sce_colnames <- lapply(sce_list, function(x) colnames(colData(x)))
  for (i in 1:length(sce_list)) {
    # There should be no difference:
    total_diffs <- length(setdiff(new_sce_colnames[[1]], new_sce_colnames[[i]]))
    if (total_diffs != 0) { 
      stop("Could not force colData column names to match across SCE objects. Cannot `cbind()` SCE objects.")
    }
  }

  # Subset SCEs to shared genes and ensure appropriate naming ------------------
  
  # First, obtain intersection among all SCE objects
  shared_genes <- Reduce(intersect, lapply(sce_list, rownames)) 
  
  # Now, loop over SCEs to subset each to the array of `shared_genes`
  #  At the same time, we also update the rowData column names to be unique across SCEs
  #  COMMENTED OUT: We also update the colData rownames to be be unique, on the off-chance we have repeated barcodes
  for (i in 1:length(sce_list)){
    
    # Subset to shared genes
    sce_list[[i]] <- sce_list[[i]][shared_genes,]
      
    # Add relevant sample IDs to rowData column names
    colnames(rowData(sce_list[[i]])) <- paste0(colnames(rowData(sce_list[[i]])), "-", sample_ids[i])
    
    # Add relevant sample IDs to colData row names
    #rownames(colData(sce_list[[i]])) <- paste0(rownames(colData(sce_list[[i]])), "-", sample_ids[i])
    
  }
  
  # Combine SCE objects with `cbind()` -----------------------------------------
  combined_sce <- do.call(cbind, sce_list)
  
  
  # Return updated sce_list and combined SCE object ----------------------------
  return(
    list(
      sce_list_updated = sce_list,
      combined_sce = combined_sce
    )
  )
}







