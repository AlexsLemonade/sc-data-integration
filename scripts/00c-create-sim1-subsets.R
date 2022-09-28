# Script used to create several subsetted simulation datasets derived from `sim1`


# Option descriptions: 
# 
# --sce_dir: Path to the folder where all SCE objects are saved locally 
# --s3_sce_bucket: Bucket on S3 where SCE objects are stored
# --copy_s3: indicates whether or not to copy existing SCE file from S3 first. 
#   To copy files use `--copy_s3`
# --overwrite: Indicates whether or not to redo sim1 subsetting and
#   overwrite any existing SCE files. To overwrite use `--overwrite`

# load the R project by finding the root directory using `here::here()`
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
    default = file.path(project_root, "sample-info", "scib-simulated-processed-libraries.tsv"),
    help = "path to metadata file listing all libraries that are to be converted"
  ),
  make_option(
    opt_str = c("--sce_dir"),
    type = "character",
    default = file.path(project_root, "data", "scib_simulated", "sce"),
    help = "path to folder where all input and output SCE objects are saved"
  ),
  make_option(
    opt_str = c("--s3_sce_bucket"),
    type = "character",
    default = "s3://sc-data-integration/scib_simulated_data/sce",
    help = "Bucket on s3 where SCE objects are stored"
  ),
  make_option(
    c("-c", "--copy_s3"),
    action = "store_true",
    help = "indicates whether or not to copy existing SCE file from S3 first. 
    To copy files use `--copy_s3`"
  ),
  make_option(
    c("-o", "--overwrite"),
    action = "store_true",
    help = "indicates whether or not to redo sim1 subsetting and 
      overwrite any existing SCE files. To overwrite use `--overwrite`"
  )
)



# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# File set up ------------------------------------------------------------------
total_sim1_batches <- 6 # there are 6 batches in sim1

# Define path to sim1 files
sim1_dir <- file.path(opt$sce_dir, "sim1")


# Ensure `sim1` files are present. There should be 6 batches.
sim1_filenames <- list.files(sim1_dir, pattern = "^sim1_", full.names=TRUE)

check_sim_files <- function(filename) {
  # Helper function for checking presence of properly-named sim SCE files
  if (!file.exists(file.path(filename))) {
    stop(glue::glue("Missing {filename}. Cannot subset `sim1` data. 
                    Make sure `00b-obtain-sim.R has been successfully run.")
    )
  }
}
purrr::walk(sim1_filenames, check_sim_files)


  
## Define cell group subsets ----------------------------------------

# The following tibbles 
# `celltypes` indicate which cell types should be retained. See discussion:
# https://github.com/AlexsLemonade/sc-data-integration/issues/151

# sim1a: All batches share all cell types (1, 2, 4). Batch 3 is removed.
sim1a_retain_celltypes <- tibble::tribble(
  ~batch, ~celltypes,
  1, c(1, 2, 4),
  2, c(1, 2, 4),
  4, c(1, 2, 4),
  5, c(1, 2, 4),
  6, c(1, 2, 4),
) 

# sim1b: Cell types are not shared across batches, but every cell type is in at
#  least 2 samples with the ability to "link" all samples
sim1b_retain_celltypes <- tibble::tribble(
  ~batch, ~celltypes,
  1, c(1, 3, 5),
  2, c(1, 2, 7),
  3, c(4, 5, 6),
  4, c(2, 3, 4),
  5, c(5, 6, 7),
  6, c(1, 2, 3),
) 

# sim1c: Cell types are not shared across batches, and some sets of samples are 
#  disjoint, without cell types that "link" all samples
sim1c_retain_celltypes <- tibble::tribble(
  ~batch, ~celltypes,
  1, c(4, 5, 6),
  2, c(4, 6, 7),
  3, c(5, 6),
  4, c(1, 2, 3),
  5, c(1, 2), 
  6, c(3, 4)
) 

# Grab existing SCE from S3 if specified -------------------------------


# Create filenames for sim1a, b, and c
sim1a_filepaths <- stringr::str_replace_all(sim1_filenames[-3], "sim1", "sim1a") # 1a does not use batch 3
sim1b_filepaths <- stringr::str_replace_all(sim1_filenames, "sim1", "sim1b")
sim1c_filepaths <- stringr::str_replace_all(sim1_filenames, "sim1", "sim1c")
all_filepaths <- c(sim1a_filenames, sim1b_filenames, sim1c_filenames)
all_filenames <- basename( all_filepaths )


# REVIEWERS!
#if(!is.null(opt$copy_s3)){
#  # copy over any existing SCEs
#  aws_includes <- paste("--include '", all_filenames, "'", sep = '', collapse = ' ')
#  sync_call <- paste('aws s3 cp', opt$s3_sce_bucket, opt$sce_dir,
#                     '--exclude "*"', aws_includes, sep = " ")
#  system(sync_call, ignore.stdout = TRUE)
#  
#}




# Create the sim1a, sim1b, and sim1c files ----------------

subset_export_sce <- function(input_sce_file, output_sce_file, celltypes_df) {
  # Helper function to read in, subset, and export a sim1 SCE file
  # input_sce_file: input file to read SCE to subset
  # output_sce_file: output file to write subsetted SCE
  # celltypes_df: data frame of which cell types to retain
  
  sce <- readr::read_rds(input_sce_file)
  
  # grab the batch index
  batch_index <- stringr::word(basename(input_sce_file), 2,sep = "_") %>%
    stringr::str_replace("Batch", "") %>%
    as.numeric()
  
  # Define array of celltypes to keep
  keep_celltypes <- paste0("Group", 
                           celltypes_df %>% 
                             dplyr::filter(batch == batch_index) %>%
                             dplyr::pull(celltypes) %>%
                             magrittr::extract2(1)
  )
  
  # Subset to only those batches
  cells_to_keep <- which(sce$celltype %in% keep_celltypes)
  sce_subsetted <- sce[, cells_to_keep]
  
  # Export sce_subsetted to RDS file
  readr::write_rds(sce_subsetted, output_sce_file)
}

# If we are not overwriting, then we only want to run files that are missing
# If we are overwriting, then we want to run all filenames
if (is.null(opt$overwrite)) {
  # Missing files only
  missing_sim1a <- sim1a_filepaths[!file.exists(sim1a_filepaths)]
  missing_sim1b <- sim1b_filepaths[!file.exists(sim1b_filepaths)]
  missing_sim1c <- sim1c_filepaths[!file.exists(sim1c_filepaths)]
} else {
  # All files
  missing_sim1a <- sim1a_filepaths
  missing_sim1b <- sim1b_filepaths
  missing_sim1c <- sim1c_filepaths
}

# Subset and export all SCE objects for each set of files that need to be created
purrr::walk2(sim1_filenames[-3], missing_sim1a, subset_export_sce, sim1a_retain_celltypes) # no batch 3 in sim1a
purrr::walk2(sim1_filenames, missing_sim1b, subset_export_sce, sim1b_retain_celltypes)
purrr::walk2(sim1_filenames, missing_sim1c, subset_export_sce, sim1c_retain_celltypes)


# REVIEWERS, help!
# Copy back to S3
#all_written_files <- c(missing_sim1a, missing_sim1b, missing_sim1c)
#aws_includes <- paste("--include '", all_written_files, "'", sep = '', collapse = ' ')
#sync_call <- paste('aws s3 sync', opt$sce_dir, opt$s3_sce_bucket, 
#                       '--exclude "*"', aws_includes, sep = " ")
#system(sync_call, ignore.stdout = TRUE)

