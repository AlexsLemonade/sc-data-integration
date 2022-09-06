
#' Create a single UMAP plot colored by a provided column in the sce object
#'
#' @param sce SCE to grab UMAP embeddings from for plotting
#' @param cell_label_column Name of column to use for coloring cells in the UMAP
#' @param umap_name Name or UMAP embeddings (e.g. "UMAP" or "fastmnn_UMAP")
#' @param plot_colors Vector of colors to use for labeling cells. Must be 
#'   equivalent to the number of categories present in the cell_label_column.
#' @param plot_title Title to use for the plot
#'
#' @return ggplot2 object containing UMAP
#'
plot_umap_panel <- function(sce,
                            cell_label_column,
                            umap_name,
                            plot_colors,
                            plot_title = NULL){
  
  # check that plot_colors is equal to the number of categories present in the cell label column 
  color_categories <- unique(colData(sce)[,cell_label_column])
  if(length(plot_colors) != length(color_categories)){
    stop("Number of colors provided must be equal to the number of categories used to classify cells in 
         the specified cell_label_column.")
  }
  
  # create umap and label with provided cell label column 
  umap <- scater::plotReducedDim(sce, 
                                 dimred = umap_name,
                                 colour_by = cell_label_column,
                                 point_size = 0.1,
                                 point_alpha = 0.4) +
    scale_color_manual(values = plot_colors) +
    # relabel legend and resize dots
    guides(color = guide_legend(title = cell_label_column,
                                override.aes = list(size = 1)))+
    ggtitle(plot_title)
  
  return(umap)
}

  
#' Create faceted UMAP comparing pre and post integration
#'
#' @param merged_sce SCE object containing all libraries merged into a single SCE object 
#'   prior to any integration
#' @param integrated_sce SCE object containing the integrated object 
#' @param integration_method The name of the method that was used for integration
#'    to create `integrated_sce`. One of: fastMNN, harmony, rpca, cca, scvi, or scanorama
#' @param group_name Name to use to describe all libraries grouped together in integrated object
#' @param cell_label_column Column to use for labeling cells in UMAP
#' @param max_celltypes Maximum number of cell types to visualize during plotting.
#'   If the `cell_label_column == celltype` and contains more than the specified `max_celltypes`
#'   then only the top N cell types, where N is `max_celltypes`, will be labeled in the UMAP and all other cells will be
#'   labeled with "other". Default is 5.
#'
#' @return Combined ggplot containing a UMAP for both the unintegrated and integrated dataset
#'   with cells colored by the specified `cell_label_column`
#'   
plot_integration_umap <- function(merged_sce,
                                  integrated_sce,
                                  integration_method,
                                  group_name,
                                  cell_label_column,
                                  max_celltypes = 5) {
  
  # check that column to label cells by is present in colData
  coldata_names <- intersect(colnames(colData(merged_sce)), colnames(colData(integrated_sce)))
  if(!cell_label_column %in% coldata_names){
    stop("Provided cell_label_column should be present in both the colData of the merged SCE and integrated SCE.")
  }
  
  # if using celltype, we only want to label the top `max_celltypes`
  if(cell_label_column == "celltype"){
    # only need to relabel if > `max_celltypes` exist
    num_celltypes <- length(unique(colData(merged_sce)[,cell_label_column]))
    if(num_celltypes > max_celltypes){
      
      merged_coldata_df <- colData(merged_sce) %>%
        as.data.frame()
    
      # select top `max_celltypes` cell types based on frequency 
      selected_celltypes <- merged_coldata_df[,cell_label_column] %>%
        table() %>%
        as.data.frame() %>%
        dplyr::arrange(desc(Freq)) %>% 
        dplyr::top_n(max_celltypes) %>%
        dplyr::pull(".") %>%
        as.character()
      
      # if not in top cell types set to "other" for both merged and integrated SCE
      merged_coldata_df <- merged_coldata_df %>%
        # first label everything outside of selected celltypes as other then if NA convert back to NA 
        dplyr::mutate(new_celltype = dplyr::if_else(celltype %in% selected_celltypes, celltype, "other"),
                      celltype = dplyr::if_else(!is.na(celltype), new_celltype, NA_character_)) %>%
        dplyr::select(-new_celltype)
      
      colData(merged_sce) <- DataFrame(merged_coldata_df)
      
      integrated_coldata_df <- colData(integrated_sce) %>%
        as.data.frame() %>%
        dplyr::mutate(new_celltype = dplyr::if_else(celltype %in% selected_celltypes, celltype, "other"),
                      celltype = dplyr::if_else(!is.na(celltype), new_celltype, NA_character_))
      
      colData(integrated_sce) <- DataFrame(integrated_coldata_df)
    
    }
  }
  
  num_colors <- length(unique(merged_sce[[cell_label_column]]))
  plot_colors <- rainbow(num_colors)
  
  pre_integration_umap <- plot_umap_panel(sce = merged_sce,
                                          cell_label_column,
                                          umap_name = "UMAP",
                                          plot_colors,
                                          plot_title = paste(group_name, "Pre-Integration"))
  
  post_integration_umap <- plot_umap_panel(sce = integrated_sce,
                                           cell_label_column,
                                           umap_name = paste0(integration_method, "_UMAP"),
                                           plot_colors,
                                           plot_title = paste(group_name, 
                                                              "Post-Integration with", 
                                                              integration_method))
  
  combined_umap <- cowplot::plot_grid(pre_integration_umap, post_integration_umap, ncol = 1)
  
  return(combined_umap)
}





#' Helper function to prepare integration method for plotting
#'
#' @param df Data frame containing `integration_method` column to prepare
#'
#' @return Data frame with additional column `integration_method_factor` that has
#'  been relabeled and reordered for plotting
prepare_integration_method <- function(df) {
  df_updated <- df %>%
    # Rename for labeling
    dplyr::mutate(integration_method_factor = dplyr::if_else(
      integration_method == "unintegrated", "Pre-Integration", "Post-integration")
    ) %>% 
    # Ensure pre-integration comes first
    dplyr::mutate(integration_method_factor = 
                    forcats::fct_relevel(
                      integration_method_factor, 
                      "Pre-Integration")
    )
  
  return(df_updated)
}




#' Function to plot batch average silhouette width (ASW) metric
#'
#' @param asw_df Data frame containing batch ASW values calculated on both 
#'  integrated and unintegrated SCEs. Expected columns are at least 
#'  `rep`, `batch_asw`, and `integration_method`
#' @param seed for sina plot reproducibility
#' @return ggplot object 
plot_batch_asw <- function(asw_df, 
                           seed = seed) {
  
  # Set seed if given
  set.seed(seed)
  
  # Set up for plotting:
  #  integration method as a factor
  #  calculate mean silhouette widths across reps
  summarized_df <- asw_df %>%
    prepare_integration_method() %>%
    # Summarize silhouette widths
    dplyr::group_by(rep, integration_method_factor) %>%
    dplyr::summarize(
      mean_batch_asw = mean(batch_asw)
    ) %>%
    dplyr::ungroup()
  
  # Find the integration method for the plot title
  integration_method <- unique(asw_df$integration_method[asw_df$integration_method != "unintegrated"])
  
  # Make the plot
  asw_plot <- ggplot2::ggplot(summarized_df) + 
    ggplot2::aes(x = integration_method_factor,
                 y = mean_batch_asw) + 
    ggplot2::geom_violin() + 
    ggforce::geom_sina(alpha = 0.6) +
    # mean +/- SE
    ggplot2::stat_summary(color = "red", size = ggplot2::rel(0.25)) +
    ggplot2::labs(
      x = "Integration status",
      y = "Batch ASW for each replicate",
      title = glue::glue("Batch ASW after integration with {integration_method}")
    ) 

  # return the plot
  return(asw_plot)
  
} 




#' Function to plot PCA regression metric
#'
#' @param pca_df Data frame containing pca regression metrics calculated on both 
#'  integrated and unintegrated SCEs. Expected columns are at least 
#'  `rep`, `pc_batch_variance`, `pc_regression_scaled`, and `integration_method`
#' @param seed for sina plot reproducibility
#' 
#' @return A cowplot plot grid of plots showing PCA regression metrics in two panels
plot_pca_regression <- function(pca_df,
                                seed = seed) {
  
  # Set seed if given
  set.seed(seed)
  
  # Set up for plotting
  pca_df <- pca_df %>%
    prepare_integration_method() 
  
  # Find the integration method for the plot title
  integration_method <- unique(pca_df$integration_method[pca_df$integration_method != "unintegrated"])
  
  # Make plot in two panels via cowplot
  #  panel 1 is pc_batch_variance, and panel 2 is pc_regression_scaled
  panel1 <- ggplot2::ggplot(pca_df) + 
    ggplot2::aes(x = integration_method_factor,
                 y = pc_batch_variance) + 
    ggplot2::geom_violin() + 
    ggforce::geom_sina(alpha = 0.6) +
    # mean +/- SE
    ggplot2::stat_summary(color = "red", size = ggplot2::rel(0.25)) +
    ggplot2::labs(
      x = "Integration status",
      y = "PC batch variance for each replicate",
      title = glue::glue("PC batch variance after integration with {integration_method}")
    ) 

  
  panel2 <- ggplot2::ggplot(pca_df) + 
    ggplot2::aes(x = integration_method_factor,
                 y = pc_regression_scaled) + 
    ggplot2::geom_violin() + 
    ggforce::geom_sina(alpha = 0.6) +
    # mean +/- SE
    ggplot2::stat_summary(color = "red", size = ggplot2::rel(0.25)) +
    ggplot2::labs(
      x = "Integration status",
      y = "PC scaled regression for each replicate",
      title = glue::glue("PC scaled regression after integration with {integration_method}")
    ) 
  
  # Combine panels into 2 rows
  pca_plots <- cowplot::plot_grid(panel1, panel2, ncol = 1)
  
  # Return plot grid
  return(pca_plots)
  
} 
