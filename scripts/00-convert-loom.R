# Script used to convert data in loom format from the HCA to SCE objects
# Grabs files from AWS if not already present locally 

# import libraries
library(magrittr)
library(optparse)
library(LoomExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-m", "--metadata_file"),
    type = "character",
    default = "../sample-info/hca-library-metadata.tsv",
    help = "path to library metadata file where each row is a library"
  ),
  make_option(
    opt_str = c("-l", "--library_file"),
    type = "character",
    default = "../sample-info/hca-processed-libraries.tsv",
    help = "path to file listing all libraries that are to be converted"
  ),
  make_option(
    opt_str = c("--loom_dir"),
    type = "character",
    default = "../data/human_cell_atlas/loom",
    help = "path to folder where all loom files should be stored"
  ),
  make_option(
    opt_str = c("--sce_output_dir"),
    type = "character",
    default = "../data/human_cell_atlas/sce",
    help = "path to folder where all output sce objects should be stored"
  ),
  make_option(
    opt_str = c("--s3_bucket"),
    type = "character",
    default = "s3://ccdl-scpca-data/human_cell_atlas_data",
    help = "Bucket on s3 where data is stored"
  ),
  make_option(
    opt_str = c("--output_metadata_file"),
    type = "character",
    default = "../sample-info/hca-downstream-analysis-metadata.tsv",
    help = "Path to write out metadata file to be used as input for downstream analysis."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

######## File set up ###########################################################

# checks that provided metadata files exist
if(!file.exists(opt$metadata_file)){
  stop("--metadata_file provided does not exist.")
}

if(!file.exists(opt$library_file)){
  stop("--library_file provided does not exist.")
}

# create directories if they don't exist 
if(!dir.exists(opt$loom_dir)){
  dir.create(opt$loom_dir, recursive = TRUE)
}

if(!dir.exists(opt$sce_output_dir)){
  dir.create(opt$sce_output_dir, recursive = TRUE)
}

# read in metadata file 
all_metadata_df <- readr::read_tsv(opt$metadata_file)

# identify datasets to be converted
# read in library file and find intersection between metadata and library file 
library_id <- readr::read_tsv(opt$library_file) %>%
  dplyr::pull(library_biomaterial_id)

# get a metadata with just libraries to be processed
# add file info for sce and loom filepaths
process_metadata_df <- all_metadata_df %>%
  dplyr::filter(library_biomaterial_id %in% library_id) %>%
  dplyr::mutate(local_loom_files = file.path(opt$loom_dir, loom_file),
                # first make filename for sce file
                local_sce_file = paste0(library_biomaterial_id, "_sce.rds"),
                # then make complete path to sce file
                output_folder = file.path(opt$sce_output_dir,
                                          tissue_group,
                                          project_name),
                local_sce_path = file.path(output_folder,
                                           local_sce_file))

# create output folders for each tissue group/project for sce results 
sce_output_folders <- unique(process_metadata_df$output_folder)
for (folder in 1:length(sce_output_folders)){
  if(!dir.exists(sce_output_folders[folder])){
    dir.create(sce_output_folders[folder], recursive = TRUE)
  }
}


######### Function for converting loom files ###################################

loom_to_sce <- function(loom_file,
                        sce_file){
 
  # import loom file as SingleCellLoomExperiment
  loom <- import(loom_file, type= "SingleCellLoomExperiment")
  
  # first grab matrix from loom file 
  # convert to sparse matrix from DelayedMatrix, sparse is required for SCE
  mtx <- loom@assays@data$matrix
  counts <- as(mtx, "dgCMatrix")
  
  # create new SCE using pieces from loom file 
  sce <- SingleCellExperiment(list(counts = counts),
                              colData = colData(loom),
                              rowData = rowData(loom),
                              metadata = metadata(loom)) 
  
  # save sce file 
  readr::write_rds(sce, sce_file)
  
  return(sce)
}

######## Grab loom and convert to SCE ##########################################

# get list of sce paths 
local_sce_files <- process_metadata_df %>%
  dplyr::pull(local_sce_path)

# list of missing sce files 
missing_sce_files <- local_sce_files[which(!file.exists(local_sce_files))]

# check if any libraries are missing a corresponding sce 
if(length(missing_sce_files) != 0){
  
  # list of all loom files that correspond to missing sce files
  loom_missing_sce <- process_metadata_df %>%
    dplyr::filter(local_sce_path %in% missing_sce_files) %>%
    dplyr::pull(local_loom_files)
  
  # if any loom files don't exist for missing SCE's then grab those from AWS S3
  if(!any(file.exists(loom_missing_sce))){
    
    # get list of folders inside s3 directory to include in copying
    aws_local_copy <- loom_missing_sce[which(!file.exists(loom_missing_sce))]
    aws_includes <- paste("--include '", aws_local_copy, sep = '', collapse = ' ')
    
    # build one sync call to copy all missing loom files 
    sync_call <- paste('aws s3 cp', opt$s3_bucket, loom_missing_sce, 
                       '--exclude "*"', includes, '--recursive', sep = " ")
    
    system(sync_call, ignore.stdout = TRUE)
  }
  
  # convert to sce objects and write files 
  sce_list <- purrr::map2(loom_missing_sce,
                          missing_sce_files, 
                          loom_to_sce)
}

###### Output metadata file ####################################################

# create downstream analysis input file
process_metadata_df %>%
  dplyr::select(sample_biomaterial_id,
                library_biomaterial_id,
                local_sce_path) %>%
  dplyr::mutate(filetering_method = "miQC") %>%
  dplyr::relocate(local_sce_path, .after = filtering_method) %>%
  readr::write_tsv(opt$output_metadata_file, col_names = FALSE)

