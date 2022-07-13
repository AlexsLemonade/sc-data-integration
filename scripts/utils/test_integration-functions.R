# Script to test integration functions
renv::load(here::here())

library(magrittr) # pipe
source(
  file.path(
    here::here(),
    "scripts", 
    "utils",
    "integration-functions.R")
)


hca_metadata <- readr::read_tsv(
  here::here("sample-info",
             "hca-processed-libraries.tsv")
  )
sce_dir <- here::here("results", 
                      "human_cell_atlas", 
                      "scpca-downstream-analyses")
tissue <- "kidney"
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







