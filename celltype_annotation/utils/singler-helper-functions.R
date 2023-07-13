# This function is used to create a stacked bar chart for all SingleR references of a given label type 
# label type is either label.fine, label.main, or label.ont

celltype_barchart <- function(celltype_assignments,
                              label_type){
  
  top_celltypes <- celltype_assignments |>
    dplyr::count(reference_factor, celltype_top) |>
    tidyr::drop_na()
  
  ggplot(top_celltypes, aes(x = reference_factor, fill = celltype_top, y = n, label = n)) + 
    geom_bar(position = "stack", stat = "identity") +
    labs(
      x = "Reference dataset",
      y = "Number of cells", 
      fill = "Cell type",
      title = label_type
    ) +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    guides(x = guide_axis(angle = 90)) 
  
}

# create a table of the celltype assignments for each reference (e.g. HumanPrimaryCellAtlasData)
celltype_table <- function(celltype_assignments){
  
  all_celltypes <- celltype_assignments |> 
    dplyr::count(reference, celltype) |> 
    dplyr::ungroup() |>
    tidyr::spread(reference, n)

  return(all_celltypes)

}


# creates a plot showing the distribution of the desired metric with reference on the x-axis
# color by whether or not labels are pruned or not
# facet by label type (main, fine, or ont)
# we will use this function for looking at delta, median delta, and score distributions
plot_metric_distribution <- function(all_results,
                                     metric_type){
  
  metric_plot <- ggplot(all_results, aes_string(x = reference, y = metric_type, group = reference, color = pruned)) +
    ggforce::geom_sina(size = 0.1, alpha = 0.2) +
    stat_summary(
      aes(group = reference),
      color = "black",
      # median and quartiles for point range
      fun = "median",
      fun.min = function(x) {
        quantile(x, 0.25)
      },
      fun.max = function(x) {
        quantile(x, 0.75)
      },
      geom = "pointrange",
      position = position_dodge(width = 0.9),
      size = 0.1
    ) +
    facet_wrap(vars(label_type))+
    guides(x = guide_axis(angle = 90),
           colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
  
  return(metric_plot)
  
}


# function to calculate the median delta for each singleR result
get_median_delta <- function(singler_result){
  
  # grab all the scores
  scores <- singler_result$scores |>
    as.matrix()
  
  # calculate the median delta (max score - median of all scores) 
  median_delta <- rowMaxs(scores) - rowMedians(scores)
  
  # add median into the singler result and return as data frame
  median_delta_result <- singler_result |>
    as.data.frame() |>
    dplyr::select(labels, delta.next, pruned.labels) |>
    dplyr::mutate(median_delta = median_delta)
  
  return(median_delta_result)
}
