library(magrittr)

# load the R project
project_root <- here::here()
renv::load(project_root)

# set up file paths 
celltype_files_directory <- file.path(project_root, "data", "human_cell_atlas", "cell_type")
project_metadata_file <- file.path(project_root, "sample-info", "hca-project-metadata.tsv")
processed_libraries_file <- file.path(project_root, "sample-info", "hca-processed-libraries.tsv")

combined_celltype_file <- file.path(project_root, "sample-info", "hca-celltype-info.tsv")

# read in project metadata
project_metadata_df <- readr::read_tsv(project_metadata_file) %>%
  tidyr::drop_na()

# append full path to cell type specific file for each project 
project_metadata_df <- project_metadata_df %>%
  dplyr::mutate(celltype_filepaths = file.path(celltype_files_directory,
                                               tissue_group,
                                               project_name,
                                               celltype_filename))

# create a tibble with sample/ library ids that have been processed 
processed_libraries_df <- readr::read_tsv(processed_libraries_file) %>%
  dplyr::select(library_biomaterial_id, sample_biomaterial_id)

# 1M immune cells --------------------------------------------------------------

# read in cell type information 
immune_cells_df <- project_metadata_df %>%
  dplyr::filter(project_name == "1M_Immune_Cells") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_csv() %>%
  dplyr::filter(cell_suspension.biomaterial_core.biomaterial_id %in% processed_libraries_df$library_biomaterial_id) %>%
  dplyr::select(library_biomaterial_id = cell_suspension.biomaterial_core.biomaterial_id,
                celltype =  annotated_cell_identity.text,
                barcode) %>%
  dplyr::mutate(project = "1M_Immune_Cells") %>%
  dplyr::left_join(processed_libraries_df)

# HumanTissueTcellActivation ---------------------------------------------------

# read in cell type information 
tcell_df <- project_metadata_df %>%
  dplyr::filter(project_name == "HumanTissueTcellActivation") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_csv() %>% 
  dplyr::filter(specimen_from_organism.biomaterial_core.biomaterial_id %in% processed_libraries_df$sample_biomaterial_id) %>%
  dplyr::select(sample_biomaterial_id = specimen_from_organism.biomaterial_core.biomaterial_id,
                celltype =  annotated_cell_identity.text,
                barcode) %>%
  dplyr::mutate(project = "HumanTissueTCellActivation") %>%
  dplyr::left_join(processed_libraries_df)

# SubstantiaNigra --------------------------------------------------------------

# Need to convert library ID from submitter ID in metadata file to SRA ID present on HCA
# "SRX7129196" = "C1B"
# "SRX7129192" = "C3"

submitter_ids <- c("C1B", "C3")

substantia_nigra_df <- project_metadata_df %>%
  dplyr::filter(project_name == "HumanBrainSubstantiaNigra") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_tsv() %>%
  dplyr::filter(Library %in% submitter_ids) %>%
  dplyr::mutate(Library = ifelse(Library == "C1B", "SRX7129196", "SRX7129192"),
                barcode = stringr::word(`Sample_Names(Library_Barcode)`, -1, sep = "_")) %>%
  dplyr::select(library_biomaterial_id = Library,
                celltype = Level_1_cell_type,
                barcode) %>%
  dplyr::mutate(project = "HumanBrainSubstantiaNigra") %>%
  dplyr::left_join(processed_libraries_df)

# Combine into one tibble ------------------------------------------------------

all_celltype_df <- dplyr::bind_rows(list(
  immune_cells_df,
  tcell_df,
  substantia_nigra_df))

readr::write_tsv(all_celltype_df, combined_celltype_file)
