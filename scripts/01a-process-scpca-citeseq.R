# Script used to process (filter and normalize) CITE-seq data (ADT counts) in the altExp slot of
# a given ScPCA sce object. This script should be run *AFTER* sce's are processed with the
# `scpca-downstream-analyses` pipeline.

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries whose SCE files _may_ need to be processed.
# --filtered_sce_dir: Path to the folder where all SCE files are stored  
# --citeseq_processed_sce_dir: Path to folder where processed SCE files will be stored
# --citeseq_name: Name of the altExp slot containing CITE-seq
# --adt_threshold: Count threshold for removing ADTs
# --repeat_processing: Flag to repeat processing if output already exists. 

#load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries and source functions
suppressPackageStartupMessages({
  library(magrittr)
  library(optparse)
  library(SingleCellExperiment)
})
source(file.path(project_root, "scripts", "utils", "integration-helpers.R"))

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--library_file"),
    type = "character",
    default = file.path(project_root, "sample-info", "scpca-processed-libraries.tsv"),
    help = "path to metadata file listing all libraries to potentially process, where the column
            `has_CITEseq` is TRUE or FALSE depending if there are CITE-seq counts to process."
  ),
  make_option(
    opt_str = c("--scpca_downstream_analyses_dir"),
    type = "character",
    default = file.path(project_root, "results", "scpca", "scpca-downstream-analyses"),
    help = "path to folder where all input sce files created by the scpca_downstream_analyses pipeline are stored"
  ),
  make_option(
    opt_str = c("--citeseq_processed_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "scpca", "citeseq_sce"),
    help = "path to folder where all output sce files are stored with normalized CITE-seq"
  ),
  make_option(
    opt_str = c("--citeseq_name"),
    type = "character",
    default = "CITEseq",
    help = "Name of the altExp slot in the sce object that holds the CITE-seq data"
  ),
  make_option(
    opt_str = c("--adt_threshold"),
    type = "integer",
    default = 1, 
    help = "Threshold of ADT counts allowed. Only ADTs above the threshold are retained;
           If an ADT's detected and mean counts are <= threshold, that ADT will be removed. 
           Default value is 1."
  ),
  make_option(
    opt_str = c("--repeat_processing"),
    default = FALSE,
    action = "store_true",
    help = "Indicates whether or not to repeat CITE-seq processing steps even if 
      output SCE objects already exist"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# File set up ------------------------------------------------------------------

# check that provided library file exists
if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# identify relevant libraries with CITE-seq data to process
library_metadata_df <- readr::read_tsv(opt$library_file)
citeseq_libraries_df <- library_metadata_df %>%
  dplyr::filter(has_CITEseq == TRUE) %>%
  dplyr::mutate(input_sce_dir = file.path(opt$scpca_downstream_analyses_dir, 
                                          sample_biomaterial_id)
  )


# Define input files
# note - input and output files are sorted to ensure they are in the same order.
input_sce_files <- sort(find_input_sce_files(citeseq_libraries_df$library_biomaterial_id, opt$scpca_downstream_analyses_dir))


# Define the output files and make sure directories exist
output_sce_files <- sort(file.path(opt$citeseq_processed_sce_dir, 
                                   citeseq_libraries_df$sample_biomaterial_id, 
                                   paste0(citeseq_libraries_df$library_biomaterial_id, "_sce.rds")
))

# quick length check:
if (!(length(input_sce_files) == length(output_sce_files))) {
  stop("Error: Could not prepare data for CITEseq processing.")
}

# helper function:
create_dir <- function(dir_name) {
  if(!dir.exists(dir_name)){
    dir.create(dir_name, recursive = TRUE)
  }
}

# Make sure each output directory exist as well
purrr::walk(dirname(output_sce_files), create_dir)

 
# Process the CITE-seq present ------------------------------------------

# Filters and normalizes ADT counts in CITE-Seq altExps
# Overall reference: http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html#normalization
process_citeseq_counts <- function(input_sce, 
                                   output_sce, 
                                   citeseq_name = opt$citeseq_name) {
  # input_sce: Path to input SCE file
  # output_sce: Path to output SCE file
  # citeseq_name: Name of CITE-seq AltExp slot
  
  #print(input_sce)
  sce <- readr::read_rds(input_sce)
  starting_cell_count <- ncol(sce) # number of starting cells, for comparing back after filtering

  # Only run if file does not exist or repeat is TRUE
  if (  !(file.exists(output_sce)) | opt$repeat_process == TRUE ) {
    
    ##### Perform QC ####
    # http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html#applying-custom-qc-filters
    retain_barcodes <- DropletUtils::cleanTagCounts(altExp(sce, citeseq_name)) %>%
      tibble::as_tibble(rownames = "cell_barcode") %>%
      dplyr::filter(discard == FALSE) %>%
      dplyr::pull(cell_barcode)
  
    # Only retain ADTs with >1 counts
    adt_df <- as.data.frame(rowData(altExp(sce))) %>%
      tibble::rownames_to_column("ADT") %>%
      # discard anything at or below given threshold
      dplyr::mutate(discard = detected <= opt$adt_threshold & mean <= opt$adt_threshold)
    
    # Vector of ADTs to retain
    retain_adts <- adt_df$ADT[adt_df$discard == FALSE]
    
    # Keep only cells that are present in retain_barcodes (note that this also filters the altExp)
    sce <- sce[,retain_barcodes]
    
    # Remove low-count ADTs as well
    altExp(sce, citeseq_name) <- altExp(sce, citeseq_name)[retain_adts, ]
    
    #### Perform normalization ####
    # http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html#cite-seq-median-norm
    
    # Calculate baseline:
    baseline <- DropletUtils::ambientProfileBimodal(altExp(sce, citeseq_name))
    
    # Calculate size factors
    size_factors <- scuttle::medianSizeFactors(altExp(sce, citeseq_name), reference = baseline)

    # Error out if any size_factors are 0
    n_zero <- sum(size_factors == 0)
    if(n_zero > 0) {
      stop(
        glue::glue(
          "Error: {n_zero} cell(s) have/has an ADT size factors of zero for {basename(input_sce)}.
          Cannot perform CITE-seq processing.
          ")
      )
    }
    
    # Print warning about number of cells removed
    percent_removed <- 100* ((starting_cell_count - ncol(sce)) / starting_cell_count)
    warning(
      glue::glue("Removed {round(percent_removed, 2)}% of cells while processing ADT counts in {basename(input_sce)}.")
    )
  
    # Print warning about ADTs removed
    discard_adts <- paste(adt_df$ADT[adt_df$discard == TRUE], collapse = ", ")
    percent_removed <- 100* ((starting_cell_count - ncol(sce)) / starting_cell_count)
    warning(
      glue::glue("The following ADTs were removed due to low counts: {discard_adts}")
    )
    
    # Finally, perform normalization with the final size factors and save back to SCE
    altExp(sce, citeseq_name) <- scater::logNormCounts(altExp(sce, citeseq_name), 
                                                       size.factors = size_factors)
    
    # Add metadata columns about filtering
    metadata(sce)$citeseq_percent_cells_removed <- percent_removed
    metadata(sce)$adts_removed <- discard_adts
    
    
    # Double check we actually did get a `logcounts` assay in there
    if (!("logcounts" %in% assayNames(altExp(sce, citeseq_name)))) {
      stop("Error in CITE-seq processing: Normalized ADT counts are missing.")
    }
    
    # Double check that the cells match in RNA and CITE, just in case
    if (!(all(colnames(sce) == colnames(altExp(sce, citeseq_name))))) {
      stop("Error in CITE-seq processing: Final RNA cell barcodes don't match ADT barcodes.")
    }
    
    # Export to file
    readr::write_rds(sce, output_sce)
  }

}


# Process data depending on repeat setting:
if (all(file.exists(output_sce_files)) & is.null(opt$repeat_processing)) {
  stop("CITE-Seq data has already been processed. To overwrite files, use the `--repeat_processing` flag.")
}
if (!(is.null(opt$repeat_processing))) {
  warning("CITE-Seq data is being re-processed and files may be overwritten.")
}
purrr::walk2(input_sce_files, 
             output_sce_files,
             process_citeseq_counts)

