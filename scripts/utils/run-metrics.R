# This script contains functions used in the report that wrap metric calculations 



#' Function to calculate and visualize LISI metrics
#'
#' @param merged_sce The list of SCE objects to analyze
#' @param batch_column The grouping column of interest
#' @param lisi_type Either cLISI (celltype) or iLISI 
#'
#' @return Data frame and ggplot object of LISI results
run_lisi <- function(sce_list, batch_column, lisi_type) {
  
  # calculate for unintegrated SCE
  lisi_unintegrated <- calculate_lisi(sce_list$unintegrated,
                                      batch_column = batch_column,
                                      unintegrated=TRUE)
  
  # calculate for integrated SCEs
  lisi_integrated_list <- integration_methods %>%
    purrr::imap(~ calculate_lisi(integrated_sce = sce_list[[.x]],
                                 batch_column = batch_column, 
                                 integration_method = .x))
  
  # create combined DF 
  lisi_df <- dplyr::bind_rows(lisi_unintegrated, lisi_integrated_list)
  
  # Visualize cLISI results
  lisi_plot <- plot_lisi(lisi_df, lisi_type = lisi_type)
  
  # Return df and plot
  return(
    list(
      lisi_df = lisi_df, 
      lisi_plot = lisi_plot
    )
  )
}

#' Function to calculate and visualize ARI metrics
#'
#' @param merged_sce The list of SCE objects to analyze
#' @param batch_column The grouping column of interest
#' @param k_range Range of k values to calculate ARI across
#' @param seed Random seed for ARI calculation
#'
#' @return Data frame and ggplot object of ARI results
run_ari <- function(sce_list, batch_column, k_range, seed) {
  
  # calculate for unintegrated SCE
  ari_unintegrated <- calculate_ari(sce_list$unintegrated,
                                    batch_column = batch_column,
                                    unintegrated=TRUE, 
                                    k_range = k_range, 
                                    seed = params$seed)
  
  # calculate for integrated SCEs
  ari_integrated_list <- integration_methods %>%
    purrr::imap(~ calculate_ari(integrated_sce = sce_list[[.x]],
                                batch_column = batch_column, 
                                integration_method = .x,
                                k_range = k_range, 
                                seed = params$seed))
  
  # create combined DF 
  ari_df <- dplyr::bind_rows(ari_unintegrated, ari_integrated_list)
  
  # Visualize ARI results with "k" as facet group
  ari_plot_k <- plot_ari(ari_df, "k")
  
  # Visualize ARI results with "integration_method" as facet group
  ari_plot_integration_method <- plot_ari(ari_df, "integration_method")
  
  # Return df and plots
  return(
    list(
      ari_df = ari_df, 
      ari_plot_k = ari_plot_k, 
      ari_plot_integration_method = ari_plot_integration_method 
    )
  )
}





#' Function to calculate and visualize ASW metrics
#'
#' @param merged_sce The list of SCE objects to analyze
#' @param batch_column The grouping column of interest
#' @param seed Random seed for silhouette width calculation
#' @param by_group Logical for whether calculations should be on a per-group basis.
#'   Default is `TRUE`.
#' @param group_label Label to include in plot. Default is "Batch"
#' 
#'
#' @return Data frame and ggplot object of ARI results
run_asw <- function(sce_list, batch_column, seed, by_group = TRUE, group_label = "Batch") {
    
  # calculate for unintegrated SCE
  asw_unintegrated <- calculate_silhouette_width(sce_list$unintegrated,
                                                 batch_column = batch_column,
                                                 unintegrated=TRUE, 
                                                 seed = seed)
  
  # calculate for integrated SCEs
  asw_integrated_list <- integration_methods %>%
    purrr::imap(~ calculate_silhouette_width(integrated_sce = sce_list[[.x]],
                                             batch_column = batch_column, 
                                             integration_method = .x,
                                             seed = seed))
  
  # create combined DF 
  asw_df <- dplyr::bind_rows(asw_unintegrated, asw_integrated_list)
  
  # Visualize results
  asw_plot <- plot_asw(asw_df, seed, by_group, group_label)
  
  # Return df and plot
  return(
    list(
      asw_df = asw_df, 
      asw_plot = asw_plot
    )
  )
}