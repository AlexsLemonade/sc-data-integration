# render singleR report for all libraries within a single project

project_root <- here::here()
renv::load(project_root)

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--scpca_project_id"),
    type = "character",
    default = "SCPCP000007",
    help = "ScPCA project ID corresponding to libraries to generate SingleR reports"
  ),
  make_option(
    opt_str = c("--library_metadata"),
    type = "character",
    default = file.path(project_root, "sample-info", "scpca-processed-libraries.tsv"),
    help = "Path to a tsv file containing project ID, library ID, and sample ID."
  ),
  make_option(
    opt_str = c("--template_rmd"),
    type = "character",
    default = file.path(project_root, "celltype_annotation", "analysis", "SingleR-non-immune-ref-comparison.Rmd"),
    help = "Path to template rmd file to render."
  ),
  make_option(
    opt_str = c("--results_dir"),
    type = "character",
    default = file.path(project_root, "celltype_annotation", "analysis", "reports"),
    help = "Path to output folder to store rendered reports."
  ),
  make_option(
    opt_str = c("--local_dir"),
    type = "character",
    default = file.path(project_root, "celltype_annotation", "data"),
    help = "path to folder where all SCE objects should be stored locally"
  ),
  make_option(
    opt_str = c("--s3_prefix"),
    type = "character",
    default = "s3://nextflow-ccdl-results/scpca/processed/results",
    help = "S3 URI prefix where SCE objects containing cell type annotations are stored"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# function to generate report --------------------------------------------------

render_template <- function(project_id,
                            sample_id,
                            library_id,
                            results_dir,
                            template_rmd,
                            local_dir,
                            s3_prefix){
  
  output_filename <- glue::glue("{library_id}_singleR_analysis.html")
  
  suppressPackageStartupMessages(rmarkdown::render(
    template_rmd,
    output_file = file.path(results_dir, output_filename),
    output_dir = results_dir,
    params = list(
      s3_data_dir = glue::glue("{s3_prefix}/{project_id}"),
      local_data_dir = local_dir,
      sample_id = sample_id,
      library_id = library_id
    ),
    envir = new.env(),
    quiet = TRUE
  )) 
  
}

# Set up and checks ------------------------------------------------------------

if(!file.exists(opt$library_metadata)){
  stop("Provided library_metadata file does not exist.")
}

results_dir <- file.path(opt$results_dir, opt$scpca_project_id)
fs::dir_create(results_dir)

library_metadata <- readr::read_tsv(opt$library_metadata)

if(!(opt$scpca_project_id %in% library_metadata$project_name)){
  stop("Provided scpca_project_id not found in metadata file.")
}

# Render reports ---------------------------------------------------------------

# filter to provided project ID and make sure column names match function arguments
library_metadata |>
  dplyr::filter(project_name == opt$scpca_project_id) |>
  dplyr::select("project_id" = project_name, 
                "sample_id" = sample_biomaterial_id, 
                "library_id" = library_biomaterial_id) |>
  # render templates for each library 
  purrr::pwalk(render_template,
               results_dir = results_dir,
               template_rmd = opt$template_rmd,
               local_dir = opt$local_dir,
               s3_prefix = opt$s3_prefix)
