# Script to test integration functions
renv::load(here::here())
util_dir <- here::here("scripts", "utils")
source(file.path(util_dir, "integration-helpers.R"))
source(file.path(util_dir, "integrate-harmony.R"))
source(file.path(util_dir, "integrate-fastMNN.R"))
source(file.path(util_dir, "integrate-seurat.R"))
source(file.path(util_dir, "calculate-LISI.R"))
source(file.path(util_dir, "calculate-ARI.R"))
source(file.path(util_dir, "calculate-silhouette-width.R"))
source(file.path(util_dir, "calculate-kBET.R"))
source(file.path(util_dir, "calculate-pca-regression.R"))


library(magrittr) # pipe
library(scpcaTools) # for transformation of sce -> seurat

hca_metadata <- readr::read_tsv(
  here::here("sample-info",
             "hca-processed-libraries.tsv")
  )
sce_dir <- here::here("results",
                      "human_cell_atlas",
                      "scpca-downstream-analyses")
tissue <- "blood"
project_metadata <- hca_metadata %>%
  dplyr::filter(tissue_group == tissue) %>%
  dplyr::mutate(sce_path = file.path(sce_dir,
                                     sample_biomaterial_id,
                                     paste0(library_biomaterial_id, "_miQC_processed_sce.rds")))

library_ids <- project_metadata$library_biomaterial_id
sce1 <- readr::read_rds(project_metadata$sce_path[1])
sce2 <- readr::read_rds(project_metadata$sce_path[2])
#sce3 <- readr::read_rds(project_metadata$sce_path[3])
#sce4 <- readr::read_rds(project_metadata$sce_path[4])

# Set up sce_list
sce_list <- list(sce1, sce2) #, sce3, sce4)
names(sce_list) <- c(library_ids[1], library_ids[2]) #, library_ids[3],  library_ids[4])

# cbind it up!
combined_sce <- combine_sce_objects(sce_list,
                                    preserve_rowdata_columns = c("Gene", "ensembl_ids", "gene_names"))

# get hvg
var_genes <- perform_hvg_selection(combined_sce,
                                   n = 5000)

# run multi batch PCA and UMAP on merged object with hvg selection
combined_sce <- perform_dim_reduction(combined_sce,
                                      var_genes = var_genes,
                                      pca_type = "multi")

# single PCA
perform_dim_reduction(combined_sce,
                      var_genes = var_genes,
                      pca_type = "single")

# should fail PCA
#perform_dim_reduction(combined_sce, var_genes = var_genes, pca_type = "pca")


# use this file for testing integration; the above combined_sce object isn't properly formatted
#  after code changes
combined_sce <- readRDS("results/human_cell_atlas/merged_sce/1M_Immune_Cells_merged_sce.rds")
# Test harmony:
integrate_harmony(combined_sce, "batch", from_pca=FALSE)
integrated_sce<- integrate_harmony(combined_sce, "batch")


# Test fastMNN:
integrated_sce <- integrate_fastMNN(combined_sce, gene_list = metadata(combined_sce)$variable_genes)


# plot UMAP pre and post integration
pre_integration <- scater::plotReducedDim(combined_sce,
                                          dimred = "UMAP",
                                          colour_by = "batch")
post_integration <- scater::plotReducedDim(integrated_sce,
                                           dimred = "harmony_UMAP",
                                           colour_by = "batch")

cowplot::plot_grid(pre_integration, post_integration, ncol = 1)

# Should fail:
# integrate_harmony(combined_sce)
# integrate_harmony(combined_sce, "not_a_column")

# scanorama integration read integrated object in
anndata_dir <- file.path(here::here(
  "results",
  "human_cell_atlas",
  "integrated_anndata"
))
scanorama_integrated_file <- file.path(anndata_dir, "1M_Immune_Cells_integrated_scanorama.h5")
scanorama_integrated_sce <- zellkonverter::readH5AD(scanorama_integrated_file, skip_assays = TRUE)

# Test Seurat Integration

# Transform sce to seurat
seurat_list <- prepare_seurat_list(combined_sce)

# Integrate data using reciprocal PCA
rpca_seurat_integrated_obj <- integrate_seurat(seurat_list, reduction_method = "rpca")

# Integrate data using canonical correlation analysis
cca_seurat_integrated_obj <- integrate_seurat(seurat_list, reduction_method = "cca")


# Test score calculation
integrated_sce <- readRDS("results/human_cell_atlas/integrated_sce/1M_Immune_Cells_integrated_harmony_sce.rds")
lisi <- calculate_lisi(integrated_sce, "batch", "harmony")

batch_ari <- calculate_ari(integrated_sce, integration_method = "harmony")


celltype_ari <- calculate_ari(integrated_sce, integration_method = "harmony", batch_column = "celltype")

batch_silhouette <- calculate_silhouette_width(integrated_sce, integration_method = "harmony")

cell_silhouette <- calculate_silhouette_width(integrated_sce, integration_method = "harmony", batch_column = "celltype")

kbet <- calculate_kbet(integrated_sce, "batch", "harmony", seed = 2022,
                       # A full kBET run takes about 3 minutes; go faster for tests with only a small n=1 k0_fraction_range
                       k0_fraction_range = 0.01) 
pca_regression <- calculate_pca_regression(integrated_sce, integration_method = "harmony", seed = 2022)


# Test score calculation on unintegrated data
# uses `combined_sce` variable; need to run through line 51 above
calculate_pca_regression(combined_sce, unintegrated=TRUE, seed = 2022)
calculate_kbet(combined_sce, "batch", unintegrated=TRUE, seed = 2022,
               # A full kBET run takes about 3 minutes; go faster for tests with only a small n=1 k0_fraction_range
               k0_fraction_range = 0.01) 
calculate_batch_asw(combined_sce, unintegrated=TRUE)
calculate_batch_ari(combined_sce, unintegrated=TRUE)
calculate_lisi(combined_sce, unintegrated=TRUE)

# Test reverse ARI 
individual_sce_dir <- here::here("data", "scib_simulated", "sce")
sce_ids <- c(
  "sim1_Batch1",
  "sim1_Batch2",
  "sim1_Batch3",
  "sim1_Batch4",
  "sim1_Batch5",
  "sim1_Batch6"
)
ind_sce_files <- glue::glue("{individual_sce_dir}/{sce_ids}_sce.rds")
individual_sce_list <- purrr::map(ind_sce_files, readr::read_rds) |> 
  purrr::set_names(sce_ids)

integrated_results_dir <- here::here("results", "scib_simulated", "integrated_sce")
integrated_sce <- file.path(integrated_results_dir, "sim1_integrated_fastmnn_sce.rds")
integrated_sce <- readr::read_rds(integrated_sce)

calculate_reverse_ari(individual_sce_list,
                      integrated_sce,
                      integration_method = "fastmnn")
