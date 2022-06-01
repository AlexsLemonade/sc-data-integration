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
    opt_str = c("--s3_loom_bucket"),
    type = "character",
    default = "s3://sc-data-integration/human_cell_atlas_data/loom",
    help = "Bucket on s3 where loom data is stored"
  ),
  make_option(
    opt_str = c("--s3_sce_bucket"),
    type = "character",
    default = "s3://sc-data-integration/human_cell_atlas_data/sce",
    help = "Bucket on s3 where loom data is stored"
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
# add file info for sce filepaths
process_metadata_df <- all_metadata_df %>%
  dplyr::filter(library_biomaterial_id %in% library_id) %>%
  # first make filename for sce file
  dplyr::mutate(local_sce_file = file.path(tissue_group,
                                           project_name,
                                           paste0(library_biomaterial_id, "_sce.rds")),
                # get path to output folder needed for creating directories 
                output_folder = file.path(opt$sce_output_dir,
                                          tissue_group,
                                          project_name),
                # then make complete path to sce file
                local_sce_path = file.path(opt$sce_output_dir,
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

####### Grab existing SCE from S3 ##############################################

# grab all library ID's that should have SCE's copied over
libraries_include <- paste("--include '", "*", library_id,"'", "*", sep = '', collapse = ' ')
sync_call <- paste('aws s3 cp', opt$s3_sce_bucket, opt$sce_output_dir, 
                   '--exclude "*"', libraries_include, '--recursive', sep = " ")
system(sync_call, ignore.stdout = TRUE)

######## Grab loom and convert to SCE ##########################################

# get list of all expected sce paths 
local_sce_paths <- process_metadata_df %>%
  dplyr::pull(local_sce_path)

# list of sce files that need to be created 
missing_sce_files <- local_sce_paths[which(!file.exists(local_sce_paths))]

# if any libraries are missing a corresponding sce, then create that sce 
if(length(missing_sce_files) != 0){
  
  # first need to see if the loom file is present to create the SCE 
  # if the loom file isn't there, grab it from AWS before converting
  
  # list of all loom files that correspond to missing sce files
  loom_missing_sce <- process_metadata_df %>%
    dplyr::filter(local_sce_path %in% missing_sce_files) %>%
    dplyr::pull(loom_file)
  
  # construct full loom path 
  loom_file_paths <- file.path(opt$loom_dir, loom_missing_sce)
  
  # if any loom files don't exist for missing SCE's then grab those from AWS S3
  if(!all(file.exists(loom_file_paths))){
    
    # get list of folders inside s3 directory to include in copying
    aws_local_copy <- loom_missing_sce[which(!file.exists(loom_file_paths))]
    aws_includes <- paste("--include '", aws_local_copy, "'", sep = '', collapse = ' ')
    
    # build one sync call to copy all missing loom files 
    sync_call <- paste('aws s3 cp', opt$s3_loom_bucket, opt$loom_dir, 
                       '--exclude "*"', aws_includes, '--recursive', sep = " ")
    
    system(sync_call, ignore.stdout = TRUE)
  }
  
  # convert to sce objects and write files 
  sce_list <- purrr::map2(loom_file_paths,
                          missing_sce_files, 
                          loom_to_sce)
  
  # sync sce output to S3 
  all_sce_files <- unique(process_metadata_df$local_sce_file)
  aws_includes <- paste("--include '", all_sce_files, "'", sep = '', collapse = ' ')
  sync_call <- paste('aws s3 sync', opt$s3_sce_bucket, opt$sce_output_dir, 
                     '--exclude "*"', aws_includes, '--recursive', sep = " ")
  system(sync_call, ignore.stdout = TRUE)
  
}
