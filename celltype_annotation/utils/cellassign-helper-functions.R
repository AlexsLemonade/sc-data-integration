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

# function to calculate the median delta for each cellassign result
# input is full table of predictions returned from running cellassign
get_median_delta_cellassign <- function(cellassign_predictions){
  
  # grab all the scores
  preds <- cellassign_predictions |>
    dplyr::select(-barcode) |>
    as.matrix()
  
  # calculate the median delta (max score - median of all scores) 
  median_delta <- rowMaxs(preds) - rowMedians(preds)
  
  # join the median delta with barcode 
  median_delta_df <- data.frame(
    barcode = cellassign_predictions$barcode,
    median_delta = median_delta
    )
  
  return(median_delta_df)
}

# create plot to look at median delta
# input is combined data frame with assigned celltype and all delta median scores
# also indicate what to color points by, default is celltype 
# input dataframe must contain `median_delta` column to color points by 
plot_median_delta <- function(celltype_results,
                              color_group = celltype){
  
  # sina plot of delta median score for cellassign
  # color_group is on the x-axis and used for calculating stats 
  delta_plot <- ggplot(celltype_results, aes(x = {{color_group}}, y = median_delta)) +
    ggforce::geom_sina(size = 0.5, alpha = 0.2) +
    stat_summary(
      aes(group = {{color_group}}),
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
    theme_bw() + 
    guides(x = guide_axis(angle = 90))
  
  return(delta_plot)
  
}

# join cell type with probability/ prediction matrix 
# input is data frame of celltype predictions, must contain the `reference`, `barcode`, and `celltype` columns
join_cellassign_results <- function(celltype_assignments,
                                    all_predictions){
  
  # split cell type assignments data frame by reference 
  combined_predictions_results <- split(celltype_assignments, celltype_assignments$reference) |>
    # join each prediction matrix with associated cell type assignments 
    purrr::map2(all_predictions,
                \(celltype, predictions){
                  predictions |>
                    dplyr::left_join(celltype, by = "barcode") |>
                    # indicate which is the assigned cell type
                    dplyr::rename(
                      "assigned_celltype" = "celltype"
                    ) |>
                    tidyr::pivot_longer(!c("barcode", "assigned_celltype", "prediction", "reference"),
                                        names_to = "celltype",
                                        values_to = "probability") |>
                    # create a column to use for coloring points by if they are associated with the cell type that gets assigned
                    dplyr::mutate(assigned_celltype = dplyr::if_else(assigned_celltype == celltype, TRUE, FALSE))
                })
  
  return(combined_predictions_results)
}

# plot probability 
# plot the distribution of probability from cell assign 
# return plot with cell type on x-axis and every prediction plotted
# points are colored by if they are assigned to the cell type on the x-axis or not
# input is dataframe with all celltype assignments and list of all predictions
plot_probability <- function(celltype_assignments,
                             all_predictions){
 
  # combine results into one singular list of data frames 
  # one for each reference 
  combined_predictions_results <- join_cellassign_results(celltype_assignments,
                                                          all_predictions)
  
  prob_plot <- combined_predictions_results |> 
    purrr::imap(\(plot_df, ref_name) {
      ggplot(plot_df, aes(x = probability, y = celltype, fill = assigned_celltype)) +
        ggridges::geom_density_ridges(alpha = 0.5) +
        ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
        theme_bw() + 
        guides(x = guide_axis(angle = 90),
               colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
        labs(title = ref_name)
    })
  
  return(prob_plot)
   
}

# plot expression of marker genes for each assigned cell type 
# input is binary reference matrix with marker genes as rows and cell types as columns
# cell type assignments with `reference`, `barcode`, and `celltype` column
# name of reference 
# sce object to use for pulling out gene expression values 
plot_marker_gene_exp <- function(ref_mtx,
                                 celltype_assignments,
                                 ref_name,
                                 annotated_sce){
  
  # create a data frame with cell type and marker gene column
  ref_genes <- ref_mtx |>
    tidyr::pivot_longer(!ensembl_id,
                        names_to = "celltype",
                        values_to = "marker_gene")|>
    # get rid of extra rows where ref genes aren't in that cell type
    dplyr::filter(marker_gene == 1) 
  
  # split by cell type to make a list of data frames
  # one per cell type
  ref_genes <- split(ref_genes, ref_genes$celltype)
  
  all_ref_gene_exp <- purrr::imap(ref_genes, 
                                     \(ref_genes_df, celltype){
                                       
                                       # filter the results to the specified cell type
                                       results_df <- celltype_assignments |>
                                         dplyr::filter(reference == ref_name,
                                                       celltype == celltype) |>
                                         # rename to avoid confusion with ref gene cell types 
                                         dplyr::rename(assigned_celltype = "celltype") |>
                                         dplyr::select(assigned_celltype, barcode)
                                       
                                       # grab the counts for all ref genes
                                       gene_exp <- annotated_sce[ref_genes_df$ensembl_id,] |> 
                                         logcounts() |>
                                         as.matrix() |>
                                         as.data.frame() |>
                                         tibble::rownames_to_column("ensembl_id") |>
                                         tidyr::pivot_longer(!ensembl_id, 
                                                             names_to = "barcode",
                                                             values_to = "gene_expression") |>
                                         dplyr::left_join(results_df, by = c("barcode")) |>
                                         dplyr::mutate(match_celltype_assignment = dplyr::if_else(assigned_celltype == celltype, 
                                                                                                  "Assigned celltype", 
                                                                                                  "Other cells")) |>
                                         # calculate mean gene expression for each gene and cell type
                                         dplyr::group_by(match_celltype_assignment, ensembl_id) |>
                                         dplyr::summarise(mean_gene_exp = mean(gene_expression))
                                       
                                     }) |>
    dplyr::bind_rows(.id = "ref_gene_celltype") |>
    # remove any cell types in the reference that are not found in our data 
    dplyr::filter(ref_gene_celltype %in% celltype_assignments$celltype)
  
  # make a faceted plot for each ref celltype 
  # plot the mean expression of each marker gene across all annotated cell types 
  gene_exp_plot <- ggplot(all_ref_gene_exp, aes(x = match_celltype_assignment, y = mean_gene_exp)) +
    geom_violin( fill = "grey") + 
    facet_wrap(vars(ref_gene_celltype), ncol = 3) +
    stat_summary(
      aes(group = match_celltype_assignment),
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
    theme_bw() +
    theme(text = element_text(size = 16)) +
    guides(x = guide_axis(angle = 90)) +
    labs(x = "")
  
  return(gene_exp_plot)
}
