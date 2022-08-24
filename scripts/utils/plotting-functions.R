
#' Create faceted UMAP comparing pre and post integration
#'
#' @param merged_sce SCE object containing all libraries merged into a single SCE object 
#'   prior to any integration
#' @param integrated_sce SCE object containing the integrated object 
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param group_name Name to use to describe all libraries grouped together in integrated object
#' @param cell_label_column Column to use for labeling cells in UMAP
#'
#' @return Combined ggplot containing a UMAP for both the unintegrated and integrated dataset
#'   with cells colored by the specified `cell_label_column`
#'   
plot_integration_umap <- function(merged_sce,
                                  integrated_sce,
                                  integration_method,
                                  group_name,
                                  cell_label_column) {
  
  # check that column to label cells by is present in colData
  coldata_names <- intersect(colnames(colData(merged_sce)), colnames(colData(integrated_sce)))
  if(!cell_label_column %in% coldata_names){
    stop("Provided cell_label_column should be present in both the colData of the merged SCE and integrated SCE.")
  }
  
  # if using celltype, we only want to label the top 5 
  if(cell_label_column == "celltype"){
    # only need to relabel if > 5 celltypes
    num_celltypes <- length(unique(colData(merged_sce)[,cell_label_column]))
    if(num_celltypes > 5){
      
      merged_coldata_df <- colData(merged_sce) %>%
        as.data.frame()
    
      # select top 5 cell types based on frequency 
      selected_celltypes <- merged_coldata_df[,cell_label_column] %>%
        table() %>%
        as.data.frame() %>%
        dplyr::arrange(desc(Freq)) %>% 
        dplyr::top_n(5) %>%
        dplyr::pull(".") %>%
        as.character()
      
      # if not in top 5 cell types set to "other" for both merged and integrated SCE
      merged_coldata_df <- merged_coldata_df %>%
        dplyr::mutate(celltype = dplyr::if_else(celltype %in% selected_celltypes, celltype, "other"))
      
      colData(merged_sce) <- DataFrame(merged_coldata_df)
      
      integrated_coldata_df <- colData(integrated_sce) %>%
        as.data.frame() %>%
        dplyr::mutate(celltype = dplyr::if_else(celltype %in% selected_celltypes, celltype, "other"))
      
      colData(integrated_sce) <- DataFrame(integrated_coldata_df)
    
    }
  }
  
  num_colors <- length(unique(merged_sce[[cell_label_column]]))
  colors <- rainbow(num_colors)
  
  pre_integration_umap <- scater::plotReducedDim(merged_sce, 
                                                 dimred = "UMAP",
                                                 colour_by = cell_label_column,
                                                 point_size = 0.1) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(title = cell_label_column,
                                override.aes = list(size = 1)))+
    ggtitle(paste(group_name, "Pre-Integration"))
  
  post_integration_umap <- scater::plotReducedDim(integrated_sce,
                                                  dimred = integrated_umap_name,
                                                  colour_by = cell_label_column,
                                                  point_size = 0.1) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(title = cell_label_column,
                                override.aes = list(size = 1)))+
    ggtitle(paste(group_name, "Post-Integration with", integration_method))
  
  combined_umap <- cowplot::plot_grid(pre_integration_umap,post_integration_umap, ncol = 1)
  
  return(combined_umap)
}
