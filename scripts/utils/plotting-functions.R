#' Setup celltypes for plotting based on the `max_celltypes` to include in celltype
#'   UMAP and ASW plots
#'
#' @param sce SCE object to process celltypes for
#' @param celltype_column Name of column in SCE containing celltypes
#' @param max_celltypes Maximum number of cell types to visualize during plotting (default 5).
#'   If there are more celltypes present than the specified `max_celltypes`, 
#'   then only the top N cell types, where N is `max_celltypes`, will be selected 
#'   and given a color
#' 
#' @return SCE object with an additional column `celltype_plot_names`, holding the
#'   reformatted celltype names to include in plots
setup_celltype_plot_names <- function(sce, 
                                      celltype_column = "celltype",
                                      max_celltypes = 5) {
  
  # how many unique cell types are there?
  num_celltypes <- length(unique(colData(sce)[,celltype_column]))
  
  # Establish new character column to hold new celltype names for use in plotting
  coldata_df <- colData(sce) %>%
    as.data.frame() %>%
    # make sure that celltype_names is a character vector and not a Factor
    # this can happen if converting from AnnData and will cause errors later on
    dplyr::mutate(celltype_plot_names = as.character(celltype))
  
  # If too many celltypes, refactor `celltype_names`
  if(num_celltypes > max_celltypes){
  
    # identify top `max_celltypes` cell types based on frequency
    selected_celltypes <- levels(forcats::fct_infreq(coldata_df[,celltype_column]))[1:max_celltypes]
    
    # if not in top cell types set to "other" for both merged and integrated SCE
    coldata_df <- coldata_df %>%
      # first label everything outside of selected celltypes as other then if NA convert back to NA
      # ensure that 'other' is last
      dplyr::mutate(celltype = as.character(celltype),
                    celltype_plot_names = dplyr::if_else(celltype %in% selected_celltypes, celltype, "other"),
                    celltype_plot_names = dplyr::if_else(is.na(celltype), NA_character_, celltype_plot_names), 
                    celltype_plot_names = forcats::fct_relevel(celltype_plot_names, "other", after = Inf))
  } 
  
  # Return coldata_df to the SCE
  colData(sce) <- DataFrame(coldata_df)
  
  # Return updated SCE
  return(sce)
}






#' Make stacked barplot of number of cells across batches from an SCE with both
#'  celltype and batch information
#'
#' @param sce_coldata colData slot of an SCE object that contains `batch_column` 
#'  and `celltype_column` columns for plotting
#' @param batch_column Name of the batch column. Default "batch"
#' @param celltype_column Name of the celltype column. Default "celltype"
#' @param plot_colors Optional vector of colors for filling barplots by cell types.
#'   If `NULL`, the plot uses the default ggplot2 palette
#' @param plot_title Title for the plot
#'
#' @return ggplot2 barplot object
plot_barplot_batch_celltype <- function(sce_coldata, 
                                        batch_column = "batch", 
                                        celltype_column = "celltype",
                                        plot_colors = NULL, 
                                        plot_title = NULL
                                        ) {
  
  batch_cell_barplot <- as.data.frame(sce_coldata) %>%
    ggplot() +
    aes_string(x = batch_column, 
               fill = celltype_column) + 
    geom_bar(color = "black", size = 0.2) +
    labs(
      x = "Batch",
      y = "Number of cells", 
      fill = "Cell type", 
      title = plot_title
    ) 
  
  # Add colors if specified
  if (!is.null(plot_colors)) {
    batch_cell_barplot <- batch_cell_barplot + 
      scale_fill_manual(values = plot_colors)
  }
  
  return(batch_cell_barplot)
}


#' Create a single UMAP plot colored by a provided column in the sce object
#'
#' @param sce SCE to grab UMAP embeddings from for plotting
#' @param cell_label_column Name of column to use for coloring cells in the UMAP
#' @param umap_name Name or UMAP embeddings (e.g. "UMAP" or "fastmnn_UMAP")
#' @param plot_colors Vector of colors to use for labeling cells. Must be
#'   equivalent to the number of categories present in the cell_label_column.
#' @param include_legend_counts Whether to include counts (like (N=###)) in legend 
#'   key labels. Default is `TRUE`.
#' @param plot_title Title to use for the plot
#' @param legend_title Legend title for colors to include in plot
#' @param seed Random seed to set prior to shuffling SCE object
#'
#' @return ggplot2 object containing UMAP
#'
plot_umap_panel <- function(sce,
                            cell_label_column,
                            umap_name,
                            include_legend_counts = TRUE,
                            plot_colors,
                            plot_title = NULL, 
                            legend_title = NULL,
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
  
  # Update `cell_label_column` to include N's, if specified
  if (include_legend_counts) {
    
    cell_label_column_sym <- rlang::ensym(cell_label_column)
    coldata_df <- as.data.frame( colData(shuffled_sce) )
    
    # Calculate counts and create a string for labeling in column `names_with_n`
    celltype_counts <- coldata_df %>%
      dplyr::count(!!cell_label_column_sym) %>%
      # new name here!
      dplyr::mutate(names_with_n = paste0(
        !!cell_label_column_sym, " (N = ", n, ")"
      )) %>%
      dplyr::select(-n)
    
    # Join into colData, and update cell_label_column string
    colData(shuffled_sce) <- coldata_df %>%
      dplyr::inner_join(celltype_counts) %>%
      DataFrame()
    cell_label_column <- "names_with_n"
  }

  # create umap and label with provided cell label column
  umap <- scater::plotReducedDim(shuffled_sce,
                                 dimred = umap_name,
                                 colour_by = cell_label_column,
                                 point_size = 0.1,
                                 point_alpha = 0.4) +
    scale_color_manual(values = plot_colors) +
    # relabel legend and resize dots
    guides(color = guide_legend(title = legend_title,
                                override.aes = list(size = 3),
                                label.theme = element_text(size = 16))) +
    theme(legend.position = "none",
          text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    labs(
      title = plot_title
    )

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
#' @param cell_label_column Column to use for labeling cells in UMAPs
#' @param include_legend_counts Whether to include counts (like (N=###)) in legend 
#'   key labels. Default is `TRUE`.
#' @param legend_title Legend title for colors to include in plot
#' @param plot_colors Optional vector of colors to use in plot.
#' @param seed Random seed to use for randomizing plotting in UMAPs
#'
#' @return Combined ggplot containing a UMAP for both the unintegrated and integrated dataset
#'   with cells colored by the specified `cell_label_column`
#'
plot_integration_umap <- function(sce,
                                  integration_method,
                                  group_name,
                                  cell_label_column,
                                  include_legend_counts = TRUE,
                                  legend_title,
                                  plot_colors = NULL,
                                  seed = NULL) {

  set.seed(seed)
  
  # check that column to label cells by is present in colData
  if(!cell_label_column %in% colnames(colData(sce))){
    stop("Provided cell_label_column should be present in the SCE object.")
  }
  
  # Define colors if not provided, or if provided check the size
  if (is.null(plot_colors)) {
    num_colors <- length(unique(sce[[cell_label_column]]))
    plot_colors <- rainbow(num_colors)    
  } else {
    if (!(length(plot_colors)) == length(unique(sce[[cell_label_column]]))) {
      stop("The number of provided colors does not match the number of labels.")
    }
  }

  if(integration_method == "unintegrated"){
    umap_name <- "UMAP"
  } else {
    # grab dim reduction name to use for plotting
    umap_name <- paste0(integration_method, "_UMAP")
  }

  umap <- plot_umap_panel(sce = sce,
                          cell_label_column,
                          umap_name = umap_name,
                          include_legend_counts, 
                          plot_colors = batch_colors,
                          plot_title = integration_method,
                          legend_title = legend_title,
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




#' Function to plot average silhouette width (ASW) metric
#'
#' @param asw_df Data frame containing silhouette width values calculated on both
#'  integrated and unintegrated SCEs. Expected columns are at least
#'  `rep`, `silhouette_width`, `silhouette_cluster`, and `integration_method`.
#' @param seed for sina plot reproducibility
#' @param by_batch Whether to take the average and color by the batch
#' @param batch_label Label to include in plot for batch, if by_batch is TRUE
#' @param plot_colors Optional vector of colors to use
#' @param label_df Optional vector frame of labels to use in the legend
#' 
#' @return ggplot object
plot_asw <- function(asw_df,
                     seed = seed, 
                     by_batch, 
                     batch_label, 
                     plot_colors = NULL, 
                     legend_labels = NULL) {

  # Set seed if given
  set.seed(seed)

  # Check that all expected columns are present in dataframe
  expected_columns <- c("integration_method", "rep", "silhouette_width", "silhouette_cluster")
  if(!all(expected_columns%in% colnames(asw_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_silhouette_width` has been run successfully.")
  }

  # Set integration method order
  asw_df_updated <- set_integration_order(asw_df)


  # Prepare and plot data by batch
  if (by_batch) {
    asw_plot <- asw_df_updated %>%
      # the `silhouette_cluster` column contains the true identity; rename for ease
      dplyr::group_by(rep, integration_method_factor, silhouette_cluster) %>%
      dplyr::summarize(
        # Use absolute value: https://github.com/AlexsLemonade/sc-data-integration/issues/149
        asw = mean(abs(silhouette_width))
      ) %>%
      dplyr::ungroup() %>%
      ggplot() +
      aes(x = integration_method_factor,
          y = asw, 
          color = silhouette_cluster) +
      ggforce::geom_sina(size = 1, alpha = 0.5,
                         position = position_dodge(width = 0.5)) +
      # add median/IQR pointrange to plot
      stat_summary(
        aes(group = silhouette_cluster),
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
      ) 
    
    if (is.null(plot_colors)) {
      asw_plot <- asw_plot + 
        ggokabeito::scale_color_okabe_ito(name = batch_label) 
    } else {
      asw_plot <- asw_plot + 
        scale_color_manual(name = batch_label, values = plot_colors, labels = legend_labels)  
    }
  } else { 
    # or without batch grouping/coloring
    asw_plot <- asw_df_updated %>%
      dplyr::group_by(rep, integration_method_factor) %>%
      dplyr::summarize(
        # Use absolute value: https://github.com/AlexsLemonade/sc-data-integration/issues/149
        asw = mean(abs(silhouette_width))
      ) %>%
      dplyr::ungroup() %>%
      ggplot() +
      aes(x = integration_method_factor,
          y = asw) +
      ggforce::geom_sina(size = 0.8, alpha = 0.7,
                         position = position_dodge(width = 0.5)) +
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
        size = 0.3
      )
    }

  # Add shared labeling
  asw_plot <- asw_plot + 
    labs(
      x = "Integration method",
      y = "Average silhouette width"
    )

  # return the plot
  return(asw_plot)

}


#' Plot adjusted rand index (ARI) across integration methods
#'
#' @param ari_df Dataframe containing the calculated ARI, either batch or celltype,
#'   after clustering across a range of values used for k (number of centers) in 
#'   k-means clustering. The dataframe must contain the following columns:
#'   "integration_method", "ari", and "k"
#' @param facet_group Group to facet plots by, can be one of "integration_method" or
#'   "k". The group not chosen for faceting will be used for the x-axis
#'
#' @return A ggplot object containing a sina plot of ARI across integration
#'   methods and values of k tested

plot_ari <- function(ari_df,
                     facet_group){

  # check that all expected columns are present in dataframe
  if(!all(c("integration_method", "ari", "k") %in% colnames(ari_df))){
    stop("Required columns are missing from input dataframe, make sure that `calculate_ari` has been run successfully.")
  }

  # check that x axes group and facet group are one of k or integration method factor
  if(!facet_group %in% c("integration_method", "k")) {
    stop("`facet_group` must be one of `integration_method` or `k`")
  }

  # set order of integration methods on axes
  ari_updated_df <- set_integration_order(ari_df)

  # define x axis, facet group variables, and plot labels
  if (facet_group == "integration_method"){
    facet_group_label <- "integration_method_factor"
    # make sure that x axes is factor if k is used
    ari_updated_df$k <- as.factor(ari_updated_df$k)
    x_axis_group <- rlang::sym("k")
    # set axes label
    x_label <- "Value of k for k-means"
  } else {
    facet_group_label <- "k"
    x_axis_group <- "integration_method_factor"
    x_label <- "Integration method"
  }


  # sina plot with batch ARI on y axis
  ari_plot <- ggplot(ari_updated_df, aes_string(x_axis_group , y = "ari")) +
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
      y = "ARI"
    ) +
    theme(axis.text.x = element_text(angle = 90))

  return(ari_plot)
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
  num_batches <- length(unique(lisi_df$batch_identity))
  
  if (lisi_type == "iLISI") {
    # normalize following the scIB method
    # https://github.com/theislab/scib/blob/067eb1aee7044f5ce0652fa363ec8deab0e9668d/scib/metrics/lisi.py#L98-L100
    lisi_df_updated <- lisi_df %>%
      dplyr::mutate(lisi_score_norm = (lisi_score-1)/(num_batches - 1))
  } else {
    # normalize following the scIB method
    # https://github.com/theislab/scib/blob/067eb1aee7044f5ce0652fa363ec8deab0e9668d/scib/metrics/lisi.py#L157-L159    
    lisi_df_updated <- lisi_df %>%
      dplyr::mutate(lisi_score_norm = (num_batches - lisi_score)/(num_batches - 1))
    
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
