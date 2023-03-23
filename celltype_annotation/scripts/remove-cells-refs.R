# Script to subset celldex references provided a list of cell types or 
# patterns found in celltype names to be removed. The cell types will be removed by matching
# the provided list with names found in the `label.main` column of the provided reference files.

project_root <- here::here()
renv::load(project_root)

library(optparse)
library(celldex)

option_list <- list(
  make_option(
    opt_str = c("--ref_names"),
    type = "character",
    default = c("BlueprintEncodeData", "HumanPrimaryCellAtlasData"),
    help = "List of celldex references to remove immune cells from."
  ),
  make_option(
    opt_str = c("--celltypes_to_remove"),
    type = "character",
    default = file.path(project_root, "celltype_annotation", "immune_celltype_patterns.txt"),
    help = "Path to a txt file containing a list of cell types or patterns present in the `label.main` column of the 
      provided references to remove from the reference. One cell type should be listed per line."
  ),
  make_option(
    opt_str = c("--ref_dir"),
    type = "character",
    default = file.path(project_root, "celltype_annotation", "references"),
    help = "path to folder where all references should be stored locally"
  ),
  make_option(
    opt_str = c("--output_prefix"),
    type = "character",
    default = "non-immune",
    help = "prefix to add to the new output files with the removed cell types."
  ),
  make_option(
    opt_str = c("--s3_prefix"),
    type = "character",
    default = "s3://scpca-references/celltype/references",
    help = "S3 URI prefix where SCE objects are stored"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# make output directory for temp ref storage
if(!dir.exists(opt$ref_dir)){
  dir.create(opt$ref_dir, recursive = TRUE)
}

# create file names from provided refnames 
ref_filenames <- glue::glue("celldex-{opt$ref_names}.rds") |>
  

# output file names 
output_filenames <- glue::glue("{opt$output_prefix}-celldex-{opt$ref_names}.rds")
output_filepaths <- file.path(opt$ref_dir, output_filenames)

# check that input list of cell types to remove has been provided 
if(!file.exists(opt$celltypes_to_remove)){
  stop("Must provide a txt file containing a list of cell types to remove.")
}

celltypes_to_remove <- readLines(opt$celltypes_to_remove)

# Sync and read in ref files ---------------------------------------------------

# sync celldex refs from aws to local directory
aws_includes <- glue::glue('--include "{ref_filenames}"') |>
  glue::glue_collapse(sep = ' ')

sync_call <- glue::glue('aws s3 sync {opt$s3_prefix} {opt$ref_dir} --exclude "*" {aws_includes}')
system(sync_call, ignore.stdout = TRUE)

# check that all files are present after syncing
local_ref_filepaths <- file.path(opt$ref_dir, ref_filenames)

missing_ref_files <- local_ref_filepaths[!(which(file.exists(local_ref_filepaths)))]
if(length(missing_ref_files)){
  stop("Provided reference names are not associated with reference files either locally or in S3 directory.")
}

# read in ref
ref_list <- local_ref_filepaths |> 
  purrr::map(readr::read_rds)

# Remove cell types from refs --------------------------------------------------

# create a regex from the list of celltypes/patterns to remove
celltype_remove_pattern <- celltypes_to_remove |> 
  paste(collapse = "|")

# subset and remove desired cell types
subset_ref_list <- purrr::walk2(ref_list, output_filepaths, 
                                \(ref, output_file){
                                  
                                  # grab index of cells to remove 
                                  removed_cell_idx <- grep(immune_celltypes, ref$label.main)
                                  
                                  # subset reference and write out to new file
                                  subset_ref <- ref[, -removed_cell_idx] |>
                                    readr::write_rds(output_file)

                                })

# Sync modified refs back to S3 ------------------------------------------------

# copy modified ref to aws 
aws_includes <- glue::glue('--include "{output_filenames}"') |>
  glue::glue_collapse(sep = ' ')

sync_call <- glue::glue('aws s3 sync {opt$ref_dir} {opt$s3_prefix} --exclude "*" {aws_includes}')
system(sync_call, ignore.stdout = TRUE)
