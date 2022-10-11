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