
#' Create a single UMAP plot colored by a provided column in the sce object
#'
#' @param sce SCE to grab UMAP embeddings from for plotting
#' @param cell_label_column Name of column to use for coloring cells in the UMAP
#' @param umap_name Name or UMAP embeddings (e.g. "UMAP" or "fastmnn_UMAP")
#' @param plot_colors Vector of colors to use for labeling cells. Must be
#'   equivalent to the number of categories present in the cell_label_column.
#' @param plot_title Title to use for the plot
#' @param seed Random seed to set prior to shuffling SCE object
#'
#' @return ggplot2 object containing UMAP
#'
plot_umap_panel <- function(sce,
                            cell_label_column,
                            umap_name,
                            plot_colors,
                            plot_title = NULL, 
                            seed = NULL){
  
  set.seed(seed)

  # check that plot_colors is equal to the number of categories present in the cell label column
  color_categories <- unique(colData(sce)[,cell_label_column])
  if(length(plot_colors) != length(color_categories)){
    stop("Number of colors provided must be equal to the number of categories used to classify cells in
         the specified cell_label_column.")
  }
  
  # randomly shuffle cells prior to plotting
  col_order <- sample(ncol(sce))
  shuffled_sce <- sce[,col_order]

  # create umap and label with provided cell label column
  umap <- scater::plotReducedDim(shuffled_sce,
                                 dimred = umap_name,
                                 colour_by = cell_label_column,
                                 point_size = 0.1,
                                 point_alpha = 0.4) +
    scale_color_manual(values = plot_colors) +
    # relabel legend and resize dots
    guides(color = guide_legend(title = cell_label_column,
                                override.aes = list(size = 3),
                                label.theme = element_text(size = 16))) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
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
#' @param seed Random seed to use for randomizing plotting in UMAPs
#'
#' @return Combined ggplot containing a UMAP for both the unintegrated and integrated dataset
#'   with cells colored by the specified `cell_label_column`
#'
plot_integration_umap <- function(sce,
                                  integration_method,
                                  group_name,
                                  cell_label_column,
                                  max_celltypes = 5, 
                                  seed = NULL) {

  set.seed(seed)
  
  # check that column to label cells by is present in colData
  if(!cell_label_column %in% colnames(colData(sce))){
    stop("Provided cell_label_column should be present in the SCE object.")
  }

  # if using celltype, we only want to label the top `max_celltypes`
  if(cell_label_column == "celltype"){
    
    # make sure that NA is actually set to NA, specifically a problem for python methods 
    colData(sce)[[cell_label_column]][which(colData(sce)[[cell_label_column]] == "NA")] <- NA_character_
    
    # only need to relabel if > `max_celltypes` exist
    num_celltypes <- length(unique(colData(sce)[,cell_label_column]))
    if(num_celltypes > max_celltypes){

      coldata_df <- colData(sce) %>%
        as.data.frame() %>%
        # make sure that celltype is a character vector and not a Factor
        # this can happen if converting from AnnData and will cause errors later on
        dplyr::mutate(celltype = as.character(celltype))

      # select top `max_celltypes` cell types based on frequency
      selected_celltypes <- coldata_df[,cell_label_column] %>%
        table() %>%
        as.data.frame() %>%
        dplyr::arrange(desc(Freq)) %>%
        dplyr::top_n(max_celltypes) %>%
        dplyr::pull(".") %>%
        as.character()

      # if not in top cell types set to "other" for both merged and integrated SCE
      coldata_df <- coldata_df %>%
        # first label everything outside of selected celltypes as other then if NA convert back to NA
        dplyr::mutate(new_celltype = dplyr::if_else(celltype %in% selected_celltypes, celltype, "other"),
                      celltype = dplyr::if_else(!is.na(celltype), new_celltype, NA_character_)) %>%
        dplyr::select(-new_celltype)

      colData(sce) <- DataFrame(coldata_df)
    }
  }

  num_colors <- length(unique(sce[[cell_label_column]]))
  plot_colors <- rainbow(num_colors)

  if(integration_method == "unintegrated"){
    umap_name <- "UMAP"
  } else {
    # grab dim reduction name to use for plotting
    umap_name <- paste0(integration_method, "_UMAP")
  }

  umap <- plot_umap_panel(sce = sce,
                          cell_label_column,
                          umap_name = umap_name,
                          plot_colors,
                          plot_title = integration_method,
                          seed = seed)

  return(umap)
}


#' Set order of integration methods
#'
#' @param metrics_df Dataframe containing desired metrics to plot, must contain
#'   a column named "integration_method"
#' @param integration_order Vector indicating the desired order of the integration methods
#'   on the axes for plotting. Default is c("unintegrated", "fastmnn", "harmony",
#'   "rpca", "cca", "scanorama", "scvi")
#'
#' @return updated dataframe with a new column, "integration_method_factor", containing the
#'   integration methods re-leveled to the desired axes_order

set_integration_order <- function(metrics_df,
                                  integration_order = c("unintegrated",
                                                        "fastmnn",
                                                        "harmony",
                                                        "rpca",
                                                        "cca",
                                                        "scanorama",
                                                        "scvi")){

  # make sure that provided dataframe contains `integration_method` as a column
  if (!"integration_method" %in% colnames(metrics_df)){
    stop("`metrics_df` is missing the column named `integration_method`.")
  }

  # check that all labels provided in the `integration_order` argument are in the integration_method column
  # if they are not present remove other integration methods from the order
  integration_methods_keep <- integration_order  %in% unique(metrics_df$integration_method)
  integration_order <- integration_order[integration_methods_keep]

  # reorder based on specified order
  updated_metrics_df <- metrics_df %>%
    dplyr::mutate(integration_method_factor = dplyr::if_else(integration_method == "unintegrated",
                                                             "Pre-Integration",
                                                             integration_method)) %>%
    dplyr::mutate(integration_method_factor = forcats::fct_relevel(integration_method, integration_order))

  return(updated_metrics_df)
}




#' Function to plot batch average silhouette width (ASW) metric
#'
#' @param asw_df Data frame containing batch ASW values calculated on both
#'  integrated and unintegrated SCEs. Expected columns are at least
#'  `rep`, `silhouette_width`, `silhouette_cluster`, `cell_name`, 
#'  and `integration_method`.
#' @param seed for sina plot reproducibility
#' @return ggplot object
plot_batch_asw <- function(asw_df,
                           seed = seed) {

  # Set seed if given
  set.seed(seed)

  # Check that all expected columns are present in dataframe
  expected_columns <- c("integration_method", "rep", "silhouette_width", "cell_name", "silhouette_cluster")
  if(!all(expected_columns%in% colnames(asw_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_batch_silhouette_width` has been run successfully.")
  }

  # Set integration method order
  asw_df_updated <- set_integration_order(asw_df)


  # Make the sina plot
  asw_plot <- asw_df_updated %>%
    # Extract library name into its own column
    dplyr::mutate(library_id = stringr::word(cell_name, -1, sep = "-")) %>%
    dplyr::group_by(rep, integration_method_factor, library_id) %>%
    dplyr::summarize(
      # Use absolute value: https://github.com/AlexsLemonade/sc-data-integration/issues/149
      mean_batch_asw = mean(abs(silhouette_width))
    ) %>%
    dplyr::ungroup() %>% 
    ggplot() +
    aes(x = integration_method_factor,
        color = library_id,
        y = mean_batch_asw) +
    ggforce::geom_sina(size = 0.8, alpha = 0.7,
                       position = position_dodge(width = 0.5)) +
    # add median/IQR pointrange to plot
    stat_summary(
      aes(group = library_id),
      color = "black",
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      position = position_dodge(width = 0.5),
      size = 0.2
    ) +
    labs(
      x = "Integration method",
      y = "Batch average silhouette width",
      color = "Batch"
    )


  # return the plot
  return(asw_plot)

}


#' Plot batch adjusted rand index (ARI) across integration methods
#'
#' @param batch_ari_df Dataframe containing the calculated batch ARI after clustering
#'   across a range of values used for k (number of centers) in k-means clustering.
#'   The dataframe must contain the following columns:
#'   "integration_method", "batch_ari", and "k"
#'@param facet_group Group to facet plots by, can be one of "integration_method" or
#'  "k". The group not chosen for faceting will be used for the x-axis
#'
#' @return A ggplot object containing a sina plot of batch ARI across integration
#'   methods and values of k tested

plot_batch_ari <- function(batch_ari_df,
                           facet_group){

  # check that all expected columns are present in dataframe
  if(!all(c("integration_method", "batch_ari", "k") %in% colnames(batch_ari_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_batch_ari` has been run successfully.")
  }

  # check that x axes group and facet group are one of k or integration method factor
  if(!facet_group %in% c("integration_method", "k")) {
    stop("`facet_group` must be one of `integration_method` or `k`")
  }

  # set order of integration methods on axes
  batch_ari_updated_df <- set_integration_order(batch_ari_df)

  # define x axis, facet group variables, and plot labels
  if (facet_group == "integration_method"){
    facet_group_label <- "integration_method_factor"
    # make sure that x axes is factor if k is used
    batch_ari_updated_df$k <- as.factor(batch_ari_updated_df$k)
    x_axis_group <- rlang::sym("k")
    # set axes label
    x_label <- "Value of k for k-means"
  } else {
    facet_group_label <- "k"
    x_axis_group <- "integration_method_factor"
    x_label <- "Integration method"
  }


  # sina plot with batch ARI on y axis
  batch_ari_plot <- ggplot(batch_ari_updated_df, aes_string(x_axis_group , y = "batch_ari")) +
    ggforce::geom_sina(size = 0.3, alpha = 0.5) +
    # facet by desired label
    facet_wrap(as.formula(paste("~", facet_group_label))) +
    # add median point to plot
    stat_summary(
      color = "red",
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      size = 0.1
    ) +
    labs(
      x = x_label,
      y = "Batch ARI"
    ) +
    theme(axis.text.x = element_text(angle = 90))

  return(batch_ari_plot)
}





#' Plot LISI score across each cell in an integrated dataset
#'
#' @param lisi_df Dataframe containing the calculated lisi (either iLISI or cLISI) scores for each cell
#'   The dataframe must contain the following columns: 
#'   "integration_method", "lisi_score", and "batch_identity" (contains either batch or cell type information)
#' @param lisi_type Either "iLISI" or "cLISI" indicating the score type
#'
#' @return A gggplot object containing both a boxplot and density plot of the given lisi score 
#'  across integration methods

plot_lisi <- function(lisi_df, lisi_type = "iLISI"){
  
  # check that all expected columns are present in dataframe 
  if(!all(c("integration_method", "lisi_score", "batch_identity") %in% colnames(lisi_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_lisi` has been run successfully.")
  }
  
  # check score type
  if (!lisi_type %in% c("cLISI", "iLISI")) {
    stop("Either `cLISI` or `iLISI` must be specified as the lisi_type.")
  }
  
  # Perform score normalization using either the number of batches or cell types depending on the score type
  if (lisi_type == "iLISI") {
    num_batches <- length(unique(lisi_df$batch_identity))
    # normalize following the scIB method
    # https://github.com/theislab/scib/blob/067eb1aee7044f5ce0652fa363ec8deab0e9668d/scib/metrics/lisi.py#L98-L100
    lisi_df_updated <- lisi_df %>%
      dplyr::mutate(lisi_score_norm = (lisi_score-1)/(num_batches - 1))
  } else {
    num_cell_types <- length(unique(lisi_df$batch_identity))
    # normalize following the scIB method
    # https://github.com/theislab/scib/blob/067eb1aee7044f5ce0652fa363ec8deab0e9668d/scib/metrics/lisi.py#L157-L159    
    lisi_df_updated <- lisi_df %>%
      dplyr::mutate(lisi_score_norm = (num_cell_types - lisi_score)/(num_cell_types - 1))
    
  }
  
  # order by median ilisi score ensuring that unintegrated will be first regardless
  lisi_df_updated <- lisi_df_updated %>%
    set_integration_order() %>%
    dplyr::mutate(integration_method_factor = forcats::fct_reorder(integration_method_factor, lisi_score, .fun = median),
                 integration_method_factor = forcats::fct_relevel(integration_method_factor, "unintegrated")
  )
  
  lisi_boxplot <- ggplot(lisi_df_updated, aes(y = lisi_score_norm, x = integration_method_factor)) +
    geom_boxplot(outlier.size = 0.25) +
    labs(
      x = "Integration method",
      y = lisi_type
    )
  
  lisi_density <- ggplot(lisi_df_updated, aes(x = lisi_score_norm, color = integration_method_factor)) +
    geom_density() +
    labs(
      color = "Integration method",
      x = lisi_type
    ) +
    ggokabeito::scale_color_okabe_ito()
  
  # create a combined boxplot and density plot
  lisi_plots <- cowplot::plot_grid(lisi_boxplot, lisi_density, ncol = 1)
  
  return(lisi_plots)
  
}




#' Plot kBET rejection rate across integration methods
#'
#' @param kbet_df Dataframe containing the calculated kBET rejection rates with the following columns:
#'   "integration_method", "kbet_stat", and "kbet_stat_type"
#'
#' @return A ggplot object containing a violin plot of kBET rejection rates across integration methods

plot_kbet <- function(kbet_df){
  
  # check that all expected columns are present in dataframe
  if(!all(c("integration_method", "kbet_stat", "kbet_stat_type") %in% colnames(kbet_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_kbet` has been run successfully.")
  }
  
  # set order of integration methods on axes
  kbet_df_updated <- set_integration_order(kbet_df)
  
  # sina plot inside violin plot with rejection rate on y axis and integration method on x axis
  ggplot(kbet_df_updated, aes(x = integration_method_factor, y = kbet_stat, color = kbet_stat_type)) +
    geom_violin(position = "dodge") +
    ggforce::geom_sina(size = 0.2, alpha = 0.5,
                       position = position_dodge(width = 0.9)) +
    # add median point to plot
    stat_summary(
      aes(group = kbet_stat_type),
      color = "black",
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      position = position_dodge(width = 0.9),
      size = 0.2
    ) +
    labs(
      x = "Integration method",
      y = "kBet rejection rate",
      color = ""
    ) +
    scale_color_discrete(labels = c("Expected", "Observed"))
  
}



#' Function to plot PCA regression metric
#'
#' @param pca_df Data frame containing pca regression metrics calculated on both
#'  integrated and unintegrated SCEs. Expected columns are at least
#'  `rep`, `pc_batch_variance`, `pc_regression_scaled`, and `integration_method`
#' @param seed for sina plot reproducibility
#'
#' @return A faceted plot showing PCA regression metrics across integration methods
plot_pca_regression <- function(pca_df,
                                seed = seed) {
  
  # Set seed if given
  set.seed(seed)
  
  
  # Check that all expected columns are present in dataframe
  expected_columns <- c("integration_method", "rep", "pc_batch_variance", "pc_regression_scaled")
  if(!all(expected_columns%in% colnames(pca_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_pca_regression()` has been run successfully.")
  }
  
  # Set up for plotting
  pca_df_updated <- pca_df %>%
    tidyr::pivot_longer(dplyr::starts_with("pc_"),
                        names_to = "metric") %>%
    dplyr::mutate(metric = ifelse(
      metric == "pc_batch_variance",
      "PC batch variance",
      "PC scaled regression"
    )) %>%
    set_integration_order()
  
  pca_reg_plot <- ggplot(pca_df_updated) +
    aes(x = integration_method_factor,
        y = value) +
    ggforce::geom_sina(size = 0.8, alpha = 0.5) +
    # facet by quantity
    facet_grid(rows = vars(metric),
               scales = "free_y") +
    # add median/IQR pointrange to plot
    stat_summary(
      color = "red",
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      size = 0.15
    ) +
    labs(
      x = "Integration method",
      y = "Metric value"
    )
  
  return(pca_reg_plot)
}
