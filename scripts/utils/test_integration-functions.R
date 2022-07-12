# Script to test integration functions and possibly serve as template for Rmd?
renv::load(here::here())

library(magrittr) # pipe
source(
  file.path(
    "scripts",
    "utils",
    "integration-functions.R")
)


hca_metadata <- readr::read_tsv(file.path("sample-info", "hca-processed-libraries.tsv"))
sce_dir <- file.path("results", "human_cell_atlas", "scpca-downstream-analyses")
project <- "KidneySingleCellAtlas"
project_metadata <- hca_metadata %>%
  dplyr::filter(project_name == project) %>%
  dplyr::mutate(sce_path = file.path(sce_dir,
                                     sample_biomaterial_id,
                                     paste0(library_biomaterial_id, "_miQC_processed_sce.rds")))

library_ids <- project_metadata$library_biomaterial_id
sce1 <- readr::read_rds(project_metadata$sce_path[1]) 
sce2 <- readr::read_rds(project_metadata$sce_path[2]) 

# Set up sce_list
sce_list <- list(sce1, sce2)
names(sce_list) <- library_ids

# cbind it up!
combined_sce <- combine_sce_objects(sce_list)







