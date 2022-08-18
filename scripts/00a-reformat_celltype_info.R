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

# Each project has individually formatted files that contain cell type information
# Each section below filters the project metadata to the specified project and grabs the 
# file path to the cell type information for that file followed by the necessary reformatting
# Each cell type file is reformatted to contain a column for 
# library_biomaterial_id, sample_biomaterial_id, barcode, celltype, and project_name

# 1M immune cells --------------------------------------------------------------

# read in cell type information 
immune_cells_df <- project_metadata_df %>%
  dplyr::filter(project_name == "1M_Immune_Cells") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_csv() %>%
  # this file only contains library Ids so filter based on that
  dplyr::filter(cell_suspension.biomaterial_core.biomaterial_id %in% processed_libraries_df$library_biomaterial_id) %>%
  dplyr::select(library_biomaterial_id = cell_suspension.biomaterial_core.biomaterial_id,
                celltype =  annotated_cell_identity.text,
                barcode) %>%
  dplyr::mutate(project = "1M_Immune_Cells") %>%
  dplyr::left_join(processed_libraries_df, by = c("library_biomaterial_id")) %>%
  dplyr::select(library_biomaterial_id,
                sample_biomaterial_id,
                project,
                barcode,
                celltype)

# HumanTissueTcellActivation ---------------------------------------------------

# read in cell type information 
tcell_df <- project_metadata_df %>%
  dplyr::filter(project_name == "HumanTissueTcellActivation") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_csv() %>% 
  # this file only contains sample ids so filter based on that
  dplyr::filter(specimen_from_organism.biomaterial_core.biomaterial_id %in% processed_libraries_df$sample_biomaterial_id) %>%
  dplyr::select(sample_biomaterial_id = specimen_from_organism.biomaterial_core.biomaterial_id,
                celltype =  annotated_cell_identity.text,
                barcode) %>%
  dplyr::mutate(project = "HumanTissueTCellActivation") %>%
  dplyr::left_join(processed_libraries_df, by = c("sample_biomaterial_id")) %>%
  dplyr::select(library_biomaterial_id,
                sample_biomaterial_id,
                project,
                barcode,
                celltype)

# SubstantiaNigra --------------------------------------------------------------

# Need to convert library ID from submitter ID in metadata file to SRA ID present on HCA
# "SRX7129196" = "C1B"
# "SRX7129192" = "C3"

submitter_ids <- c("C1B", "C3")

substantia_nigra_df <- project_metadata_df %>%
  dplyr::filter(project_name == "HumanBrainSubstantiaNigra") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_tsv() %>%
  # this file only contains the submitter specific ids that are not found on HCA
  dplyr::filter(Library %in% submitter_ids) %>%
  dplyr::mutate(Library = dplyr::case_when(Library == "C1B" ~ "SRX7129196",
                                           Library == "C3" ~ "SRX7129192"),
                # reformat barcode column originally provided
                barcode = stringr::word(`Sample_Names(Library_Barcode)`, -1, sep = "_")) %>%
  dplyr::select(library_biomaterial_id = Library,
                celltype = Level_1_cell_type,
                barcode) %>%
  dplyr::mutate(project = "HumanBrainSubstantiaNigra") %>%
  dplyr::left_join(processed_libraries_df, by = c("library_biomaterial_id")) %>%
  dplyr::select(library_biomaterial_id,
                sample_biomaterial_id,
                project,
                barcode,
                celltype)

# Oligodendrocyte_MS -----------------------------------------------------------

submitter_ids <- c("CO28", "CO14")

oligo_df <- project_metadata_df %>%
  dplyr::filter(project_name == "Oligodendrocyte_MS") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_tsv() %>% 
  # this file only contains the submitter specific ids that are not found on HCA
  dplyr::filter(Sample %in% submitter_ids) %>%
  dplyr::mutate(Sample = dplyr::case_when(Sample == "CO28" ~ "Control_CO28", 
                                          Sample == "CO14" ~ "Control_CO14"),
                # extract barcode and remove extra x at the end 
                barcode = stringr::str_remove(stringr::word(Detected, -1, sep = ":"), "x")) %>%
  dplyr::select(sample_biomaterial_id = Sample,
                celltype =  Celltypes,
                barcode) %>%
  dplyr::mutate(project = "Oligodendrocyte_MS") %>%
  dplyr::left_join(processed_libraries_df, by = c("sample_biomaterial_id")) %>%
  dplyr::select(library_biomaterial_id,
                sample_biomaterial_id,
                project,
                barcode,
                celltype)

# MultipleSclerosisLineageDiversity --------------------------------------------

submitter_ids <- c("C2", "C9")

ms_lineage_df <- project_metadata_df %>%
  dplyr::filter(project_name == "MultipleSclerosisLineageDiversity") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_tsv() %>% 
  dplyr::filter(sample %in% submitter_ids) %>%
  dplyr::mutate(sample = dplyr::case_when(sample == "C2" ~ "SRX5897026", 
                                          sample == "C9" ~ "SRX5897034"),
                # extract barcode and remove extra x at the end 
                barcode = stringr::word(cell, 1 , sep = "-"))%>%
  dplyr::select(library_biomaterial_id = sample,
                celltype =  cell_type,
                barcode) %>%
  dplyr::mutate(project = "MultipleSclerosisLineageDiversity") %>%
  dplyr::left_join(processed_libraries_df, by = c("library_biomaterial_id")) %>%
  dplyr::select(library_biomaterial_id,
                sample_biomaterial_id,
                project,
                barcode,
                celltype)

# FetalLiverHaematopoiesis -----------------------------------------------------

# submitter IDs here are slightly different "F6_kidney_CD45+", "F3_kidney_CD45+"
submitter_ids <- c("F6_kidney_CD45+", "F3_kidney_CD45+")

fetal_liver_df <- project_metadata_df %>%
  dplyr::filter(project_name == "FetalLiverHaematopoiesis") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_csv() %>%
  # extract submitter id from cell barcode 
  dplyr::mutate(submitter_id = stringr::word(`cell barcode`,1, 3, sep = "_"),
                barcode = stringr::word(`cell barcode`, -1, sep = "_")) %>%
  dplyr::filter(submitter_id %in% submitter_ids) %>%
  dplyr::mutate(library_biomaterial_id = dplyr::case_when(submitter_id == "F6_kidney_CD45+" ~ "F06-KID-3p-CD45pos",
                                                          submitter_id == "F3_kidney_CD45+" ~ "F03-KID-3p-CD45pos")) %>%
  dplyr::select(library_biomaterial_id,
                celltype = `cell labels`,
                barcode) %>%
  dplyr::mutate(project = "FetalLiverHaematopoiesis") %>%
  dplyr::left_join(processed_libraries_df, by = c("library_biomaterial_id")) %>%
  dplyr::select(library_biomaterial_id,
                sample_biomaterial_id,
                project,
                barcode,
                celltype)

# KidneySingleCellAtlas --------------------------------------------------------

kidney_df <- project_metadata_df %>%
  dplyr::filter(project_name == "KidneySingleCellAtlas") %>%
  dplyr::pull(celltype_filepaths) %>%
  readr::read_csv() %>%
  # this file only contains library Ids so filter based on that
  dplyr::filter(cell_suspension.biomaterial_core.biomaterial_id %in% processed_libraries_df$library_biomaterial_id) %>%
  dplyr::select(library_biomaterial_id = cell_suspension.biomaterial_core.biomaterial_id,
                celltype =  annotated_cell_identity.text,
                barcode) %>%
  dplyr::mutate(project = "KidneySingleCellAtlas") %>%
  dplyr::left_join(processed_libraries_df, by = c("library_biomaterial_id")) %>%
  dplyr::select(library_biomaterial_id,
                sample_biomaterial_id,
                project,
                barcode,
                celltype)

# Combine into one tibble ------------------------------------------------------

# combine all individually reformatted celltype tibbles into one tibble before saving file

all_celltype_df <- dplyr::bind_rows(
  immune_cells_df,
  tcell_df,
  substantia_nigra_df,
  oligo_df,
  ms_lineage_df,
  fetal_liver_df,
  kidney_df 
)

readr::write_tsv(all_celltype_df, combined_celltype_file)
