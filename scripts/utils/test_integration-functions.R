# Script to test integration functions
renv::load(here::here())
util_dir <- here::here("scripts", "utils")
source(file.path(util_dir, "integration-helpers.R"))
source(file.path(util_dir, "integrate-harmony.R"))
source(file.path(util_dir, "integrate-fastMNN.R"))
source(file.path(util_dir, "integrate-seurat.R"))
source(file.path(util_dir, "calculate-iLISI.R"))
source(file.path(util_dir, "calculate-batch-ARI.R"))
source(file.path(util_dir, "calculate-batch-ASW.R"))
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


# Test integration:
harmony_sce <- integrate_harmony(combined_sce, "batch", from_pca=FALSE) # counts logcounts harmony_corrected 
fastMNN_sce <-integrate_fastMNN(combined_sce) # counts logcounts fastmnn_corrected
cca_sce <- integrate_seurat(prepare_seurat_list(combined_sce), reduction_method = "cca")  # RNA SCT
rpca_sce <- integrate_seurat(prepare_seurat_list(combined_sce), reduction_method = "rpca") # RNA SCT




integrated_sce_obj <- as.SingleCellExperiment(cca_sce, 
                                              assay = "RNA") #logcounts, which comes from...


# Convert reducedDims names back to <lowercase>_<UPPERCASE> because 
#   `as.SingleCellExperiment` makes them all uppercase
reducedDimNames(integrated_sce_obj) <- stringr::str_replace_all(
  reducedDimNames(integrated_sce_obj),
  toupper(opt$seurat_reduction_method),
  opt$seurat_reduction_method
)

# Restore the original `counts` assay into integrated_sce_obj because
#  `as.SingleCellExperiment` only keeps the `logcounts` assay
counts(integrated_sce_obj) <- counts(combined_sce)







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


# Test score calculation
integrated_sce <- readRDS("results/human_cell_atlas/integrated_sce/1M_Immune_Cells_integrated_harmony_sce.rds")
lisi <- calculate_ilisi(integrated_sce, "batch", "harmony")
batch_ari <- calculate_batch_ari(integrated_sce, integration_method = "harmony")
batch_asw <- calculate_batch_asw(integrated_sce, integration_method = "harmony")
kbet <- calculate_kbet(integrated_sce, "batch", "harmony", seed = 2022,
                       # A full kBET run takes about 3 minutes; go faster for tests with only a small n=1 k0_fraction_range
                       k0_fraction_range = 0.01) 
pca_regression <- calculate_pca_regression(integrated_sce, integration_method = "harmony", seed = 2022)
