# Script to test integration functions
renv::load(here::here())
util_dir <- here::here("scripts", "utils")
source(file.path(util_dir, "integration-helpers.R"))
source(file.path(util_dir, "integrate-harmony.R"))
source(file.path(util_dir, "integrate-fastMNN.R"))
source(file.path(util_dir, "integrate-seurat.R"))

library(magrittr) # pipe
library(scpcaTools) # for transformation of sce -> seurat

hca_metadata <- readr::read_tsv(
  here::here("sample-info",
             "hca-processed-libraries.tsv")
  )
sce_dir <- here::here("results", 
                      "human_cell_atlas", 
                      "scpca-downstream-analyses")
tissue <- "brain"
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
                                    c("Gene", "ensembl_ids", "gene_names"))

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

# Test harmony:
integrate_harmony(combined_sce, "batch", from_pca=FALSE)
integrated_object <- integrate_harmony(combined_sce, "batch")

# add UMAP from harmony_PCA 
integrated_object <- perform_dim_reduction(integrated_object, 
                                           var_genes = var_genes,
                                           prefix = "harmony")

# Test fastMNN:
integrate_fastMNN(combined_sce)
# plot UMAP pre and post integration 
pre_integration <- scater::plotReducedDim(combined_sce, 
                                          dimred = "UMAP", 
                                          colour_by = "batch")
post_integration <- scater::plotReducedDim(integrated_object, 
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
  "anndata",
  "integrated_scanorama_objects"
))
scanorama_integrated_file <- file.path(anndata_dir, "1M_Immune_Cells_scanorama_integrated.h5")
scanorama_integrated_sce <- zellkonverter::readH5AD(scanorama_integrated_file,
                                                    X_name = "counts")

# Test Seurat Integration

# Transform sce to seurat
seurat1 <- scpcaTools::sce_to_seurat(sce1)
seurat2 <- scpcaTools::sce_to_seurat(sce2)

# Use `Seurat::SCTransform` for normalization
seurat1 <- Seurat::SCTransform(seurat1)
seurat2 <- Seurat::SCTransform(seurat2)

# Add batch column to the metadata
seurat1 <- Seurat::AddMetaData(seurat1, library_ids[1], "batch")
seurat2 <- Seurat::AddMetaData(seurat2, library_ids[2], "batch")

# Define seurat list of objects
seurat_list <- c(seurat1, seurat2)

# Integrate data using reciprocal PCA
rpca_seurat_integrated_obj <- integrate_seurat(seurat_list, reduction_method = "rpca", integration_dims = c(1:50))

# Integrate data using canonical correlation analysis
cca_seurat_integrated_obj <- integrate_seurat(seurat_list, reduction_method = "cca", integration_dims = c(1:50))
