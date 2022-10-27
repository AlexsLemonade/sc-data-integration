# Script used to process (filter and normalize) CITE-seq data (ADT counts) in the altExp slot of
# a given ScPCA sce object. This script should be run *AFTER* sce's are processed with the
# `scpca-downstream-analyses` pipeline.

# Option descriptions: 
# 
# --library_file: The path to the file listing all libraries whose SCE files _may_ need to be processed.
# --filtered_sce_dir: Path to the folder where all SCE files are stored  
# --citeseq_processed_sce_dir: Path to folder where processed SCE files will be stored
# --citeseq_name: Name of the altExp slot containing CITE-seq
# --repeat_processing: Flag to repeat processing if output already exists. 

#load the R project by finding the root directory using `here::here()`
project_root <- here::here()
renv::load(project_root)

# import libraries
suppressPackageStartupMessages({
  library(magrittr)
  library(optparse)
  library(SingleCellExperiment)
})

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
    default = file.path(project_root, "data", "scpca", "scpca-sownstream-analyses"),
    help = "path to folder where all input sce files created by the scpca_downstream_analyses pipeline are stored"
  ),
  make_option(
    opt_str = c("--citeseq_processed_sce_dir"),
    type = "character",
    default = file.path(project_root, "results", "scpca", "scpca-downstream-analyses_citeseq"),
    help = "path to folder where all output sce files are stored with normalized CITE-seq"
  ),
  make_option(
    opt_str = c("--citeseq_name"),
    type = "character",
    default = "CITEseq",
    help = "Name of the altExp slot in the sce object that holds the CITE-seq data"
  ),
  make_option(
    opt_str = c("--repeat_processing"),
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

# identify relevant libraries
library_metadata_df <- readr::read_tsv(opt$library_file)

# find SCE files that match library ID
library_search <- paste(opt$library_file$, collapse = "|")
all_library_files <- list.files(opt$sce_dir,
                                pattern = library_search,
                                recursive = TRUE,
                                full.names = TRUE)
# just include RDS files, otherwise HTML files will also be found
sce_files <- all_library_files[grep(pattern = ".rds", all_library_files, ignore.case = TRUE)]

# if the number of sce files is different then the library ID's find the missing files
if(length(sce_files) < length(library_ids)){
  
  libraries_found <- stringr::str_extract(sce_files, library_search)
  missing_libraries <- setdiff(library_ids, libraries_found)
  
  stop(
    glue::glue(
      "\nMissing SCE object for {missing_libraries}.
      Make sure that you have run `01-run-downstream-analyses.sh`."
    )
  )
}







# construct path to filtered RDS files and ensure they are present
# if missing, throw error to run downstream analyses 
filtered_sce_files <- library_metadata_df %>%
  
  # SCPCL000208_miQC_processed_sce.rds 
  dplyr::mutate(sce_filename = paste0(library_biomaterial_id, "_filtered.rds"),
                downstream_file_path = file.path(scpca_project_id,
                                               sample_biomaterial_id,
                                               sce_filename)) %>%
  dplyr::pull(filtered_file_path)
filtered_sce_files_input <- file.path(opt$filtered_sce_dir,
                                      filtered_sce_files)
missing_sce_files <- filtered_sce_files_input[which(!file.exists(filtered_sce_files_input))]
if(length(missing_sce_files) > 0){
  stop("Input SCE files are not present. Make sure `00d-obtain-scpca-sce.R` has been run.")
}

# Setup output:

# helper function:
create_dir <- function(dir_name) {
  if(!dir.exists(dir_name)){
    dir.create(dir_name, recursive = TRUE)
  }
}

# make sure the output directory exists
create_dir(opt$citeseq_processed_sce_dir)


# Define the output file names, using same nested directory structure as in opt$filtered_sce_dir 
filtered_sce_files_output <- file.path(opt$citeseq_processed_sce_dir,
                                       filtered_sce_files)

# Make sure each subdirectory exists as well
purrr::walk(dirname(filtered_sce_files_output), create_dir)

 
# Process the CITE-seq, where present ------------------------------------------

# For files with CITE-seq, this function processes and exports them to the output directory
# For files without CITE-seq, this function effectively copies them to the output directory 
#  with no futher processing
process_citeseq_counts <- function(input_sce, 
                                   output_sce, 
                                   citeseq_name = opt$citeseq_name) {
  # input_sce: Path to input SCE file
  # output_sce: Path to output SCE file
  # citeseq_name: Name of CITE-seq AltExp slot
  
  #print(input_sce)
  sce <- readr::read_rds(input_sce)
  starting_cell_count <- ncol(sce) # number of starting cells, for comparing back after filtering

  # Check if CITE-Seq
  if (citeseq_name %in% altExpNames(sce)) {
    # Overall reference: http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html#normalization
    
    ##### Perform QC ####
    # http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html#applying-custom-qc-filters
    retain_barcodes <- DropletUtils::cleanTagCounts(altExp(sce, citeseq_name)) %>%
      tibble::as_tibble(rownames = "cell_barcode") %>%
      dplyr::filter(discard == FALSE) %>%
      dplyr::pull(cell_barcode)
      
    # Keep only cells that are present in retain_barcodes
    # note that this also filters the altExp
    sce <- sce[,retain_barcodes]
    
    #### Perform normalization ####
    # http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html#cite-seq-median-norm
    
    # Calculate baseline:
    baseline <- DropletUtils::ambientProfileBimodal(altExp(sce, citeseq_name))
    
    # Calculate size factors
    size_factors <- scuttle::medianSizeFactors(altExp(sce, citeseq_name), reference = baseline)
    
    # Ensure that no size_factors are 0
    while(sum(size_factors == 0)) {
      
      # Remove those cells that are still 0, and start again
      sce <- sce[,size_factors != 0]
      
      baseline <- DropletUtils::ambientProfileBimodal(altExp(sce, citeseq_name))
      size_factors <- scuttle::medianSizeFactors(altExp(sce, citeseq_name), reference = baseline)
    }
    
    # Print warning about number of cells removed
    percent_removed <- 100* ((starting_cell_count - ncol(sce)) / starting_cell_count)
    warning(
      glue::glue("{round(percent_removed, 2)}% of cells while processing CITE-seq counts in {basename(input_sce)}.")
    )
  
    # Finally, perform normalization with the final size factors and save back to SCE
    altExp(sce, citeseq_name) <- scater::logNormCounts(altExp(sce, citeseq_name), 
                                                       size.factors = size_factors)
    
    # Double check we actually did get a `logcounts` assay in there
    if (!("logcounts" %in% assayNames(altExp(sce, citeseq_name)))) {
      stop("Error in CITE-seq processing: Normalized counts are missing.")
    }
    
    # Double check that the number of cells matches in RNA and CITE, just in case
    if (!(all(colnames(sce) == colnames(altExp(sce, citeseq_name))))) {
      stop("Error in CITE-seq processing: Final RNA cell barcodes don't match ADT barcodes.")
    }
    
  }
  
  # Export to file
  readr::write_rds(sce, output_sce)
}


# Process data depending on repeat setting:
if (all(file.exists(filtered_sce_files_output))) {
  if (!(is.null(opt$repeat_processing))) {
    warning("CITE-Seq data is being re-processed and files will be overwritten.")
    purrr::walk2(filtered_sce_files_input, 
                 filtered_sce_files_output,
                 process_citeseq_counts)
  } else {
    warning("CITE-Seq data has already been processed. To overwrite files, use the `--repeat_processing` flag.")
  } 
} else {
  purrr::walk2(filtered_sce_files_input, 
               filtered_sce_files_output,
               process_citeseq_counts)
}

