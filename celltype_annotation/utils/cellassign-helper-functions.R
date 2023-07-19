# create binary gene by cell matrix given a set of marker genes 

build_binary_mtx <- function(marker_genes_df, 
                             celltype_column,
                             gene_id_column,
                             rowdata_df){
 
  # create a binary matrix of genes by cell type markers 
  binary_mtx <- marker_genes_df |>
    dplyr::select(celltype_column, gene_id_column) |> 
    tidyr::pivot_wider(id_cols = gene_id_column,
                       names_from = celltype_column,
                       values_from = celltype_column,
                       values_fn = length,
                       values_fill = 0) |> 
    tibble::column_to_rownames(gene_id_column) |>
    # add a column with no marker genes 
    # cell assign will assign cells to "other" when no other cell types are appropriate 
    dplyr::mutate(other = 0)
  
  # replace length with 1
  binary_mtx[binary_mtx > 1] <- 1
  
  # merge reference with rowdata and replace gene symbols with ensembl
  binary_mtx <- binary_mtx |>
    tibble::rownames_to_column("gene_symbol") |>
    dplyr::left_join(rowdata_df) |>
    # drop any rows where there is no ensembl id
    tidyr::drop_na(ensembl_id) |>
    dplyr::select(ensembl_id, everything(), -gene_symbol)
  
  return(binary_mtx)
}

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
                                 cell_assign,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 title = NULL){
  if(is.null(title)){
    title <- ""
  }
  
  # build a matrix with reference as the rows and the celltype as the columns
  label_mtx <- table(original_assignment,
                     cell_assign,
                     useNA = "ifany") |> 
    log1p() # log transform for visual help
  
  ComplexHeatmap::pheatmap(label_mtx,
                           cluster_rows = cluster_rows, 
                           cluster_cols = cluster_cols,
                           fontsize_col = 8,
                           heatmap_legend_param = list(title = "Log(Number of cells)"),
                           main = title)
  
}
