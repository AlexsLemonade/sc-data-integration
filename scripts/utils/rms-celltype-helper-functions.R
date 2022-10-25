#' function that checks that integrated files exist, and grabs coldata  
#' Specifically for integrated Seurat objects from RMS project 
#'
#' @param integrated_file Path to RDS file containing an integrated Seurat object 
#'
#' @return Data frame containing the cell metadata and the above mentioned columns 

grab_integrated_metadata <- function(integrated_file){
  
  if(!file.exists(integrated_file)){
    stop(paste("--rms_seurat_dir must contain integrated seurat object file:", basename(integrated_file)))
  }
  
  # read in integrated object 
  integrated_seurat_obj <- readr::read_rds(integrated_file)
  
  # check that all columns are included as columns in the meta.data
  if(!metadata_cols %in% colnames(integrated_seurat_obj@meta.data)){
    stop("metadata_cols must be columns in the meta.data of the seurat object.")
  }
  
  # grab coldata as dataframe 
  coldata_df <- as.data.frame(integrated_seurat_obj@meta.data) %>%
    tibble::rownames_to_column("unique_cell_id") %>%
    dplyr::select(unique_cell_id, # sample id + barcode
                  orig.ident, # original sample id
                  cell.cluster.ids, # tumor or malignant
                  SJID, # corresponding submitter ID
                  annot.clusters) # annotated tumor cell types))
  
  # remove integrated object to save space 
  rm(integrated_seurat_obj)
  
  return(coldata_df)
  
}

#' function to read in sce object and seurat object and grab cell type information 
#' specifically for grabbing cell type information from seurat objects for RMS dataset
#' exports modified sce object 
#'
#' @param sce_processed_filepath Full path to RDS file containing normalized SCE
#' @param seurat_filepath Full path to RDS file containing individual Seurat object with 
#'   `cell.cluster.ids` column in cell metadata
#' @param library_biomaterial_id Unique ID corresponding to library 
#' @param submitter_id Unique ID associated with submitter_id from ScPCA
#' @param all_integrated_coldata Data frame from integrated seurat objects 
#'
#' @return Exports SCE object with additional celltype column and returns modified colData

add_celltype <- function(sce_processed_filepath, 
                         seurat_filepath, 
                         library_biomaterial_id, 
                         submitter_id, 
                         all_integrated_coldata){
  
  # find seurat obj 
  if(!file.exists(seurat_filepath)){
    seurat_coldata <- NULL
    glue::glue("
               No seurat object found for library ID: {library_biomaterial_id} with submitter ID: {submitter_id}.
               ")
  } else {
    # read in seurat object 
    seurat_obj <- readr::read_rds(seurat_filepath)
    # extract col data and add cell barcode column 
    seurat_coldata <- as.data.frame(seurat_obj@meta.data) %>%
      tibble::rownames_to_column("unique_cell_id") %>%
      dplyr::select(unique_cell_id,
                    "celltype" = cell.cluster.ids) %>%
      dplyr::mutate(barcode = stringr::word(unique_cell_id, -1, sep = "-"))
    
    # save space by removing seurat object 
    rm(seurat_obj)
  }
  
  # check if submitter_id in SJID column of all_coldata
  if(submitter_id %in% unique(all_integrated_coldata$SJID)){
    has_tumor_subtypes <- TRUE
  } else {
    # if not in integrated data then no subtypes for tumor cells can be found
    has_tumor_subtypes <- FALSE
  }
  
  # read in SCE object 
  sce <- readr::read_rds(sce_processed_filepath)
  sce_coldata <- as.data.frame(colData(sce)) %>%
    tibble::rownames_to_column("barcode")
  
  if(!is.null(seurat_coldata)){
    
    # join with integrated coldata only if tumor subtypes exist 
    if(has_tumor_subtypes){
      seurat_coldata <- seurat_coldata %>%
        dplyr::left_join(all_integrated_coldata) %>%
        dplyr::mutate(celltype = ifelse(celltype == "Tumor", paste(celltype, annot.clusters, sep = "_"), cell.cluster.ids)) %>%
        dplyr::select(barcode, celltype)
    }
    
    
    # join celltype information from seurat object with coldata from sce object 
    sce_coldata <- sce_coldata %>%
      dplyr::left_join(seurat_coldata)
    
  } else {
    # set celltype to NA if no seurat object is found for this submitter ID and therefore no celltype 
    sce_coldata <- sce_coldata %>%
      dplyr::mutate(celltype = NA)
  }
  # add coldata back to SCE object 
  colData(sce) <- DataFrame(sce_coldata, row.names = sce_coldata$barcode)
  
  # export rds file with annotated celltype 
  celltype_sce_filepath <- file.path(opt$celltype_sce_dir, paste0(library_biomaterial_id, "_processed_celltype.rds"))
  readr::write_rds(sce, celltype_sce_filepath)
  
  return(sce_coldata)
  
}