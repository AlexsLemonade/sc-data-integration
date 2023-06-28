
# define a function for obtaining cell type assignments from cell assign predictions
get_celltype_assignments <- function(predictions){
  
  # get a dataframe of assignment for each barcode and prediction score
  celltype_assignments <- predictions |>
    tidyr::pivot_longer(!barcode,
                        names_to = "celltype",
                        values_to = "prediction") |>
    dplyr::group_by(barcode) |>
    dplyr::slice_max(prediction, n = 1) |>
    dplyr::ungroup() 
  
  return(celltype_assignments)
}


# function for creating comparison heatmaps between two different annotations
compare_refs_heatmap <- function(original_assignment,
                                 cell_assign){
  # build a matrix with reference as the rows and the celltype as the columns
  label_mtx <- table(original_assignment,
                     cell_assign,
                     useNA = "ifany") |> 
    log1p() # log transform for visual help
  
  pheatmap::pheatmap(label_mtx,
                     cluster_rows = TRUE, 
                     width = 10,
                     fontsize_col = 8)
  
}