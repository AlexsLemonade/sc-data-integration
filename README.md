# sc-data-integration

This repo contains files, scripts, and analysis related to exploring integration of single-cell and single-nuclei data. 

## Managing R packages with `renv`

Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html). 
For `renv` to work as intended, you'll need to work within the `sc-data-integration.Rproj` project in RStudio.
You may need to run `renv::restore()` upon opening the project to ensure the `renv.lock` file is synced with the project library. 

Each time you install or use new packages, you will want to run `renv::snapshot()` to update the `renv.lock` file with any added package and dependencies necessary to run the analyses and scripts in this repo.

If there are dependencies you want to include that are not captured automatically by `renv::snapshot()`, add them to `components/dependencies.R` with a call to `library()` and an explanatory comment.
   For example, if `dplyr` were recommended but not required by a package and you wanted to make sure to include it in the lockfile, you would add `library(dplyr)` to `components/dependencies.R`.
   Then rerun `renv::snapshot()`.

## Metadata 

Inside the `sample-info` folder is metadata related to datasets used for testing data integration.

For exploring data integration, we are using test datasets that have been obtained from the [Human Cell Atlas (HCA) Data Portal](https://data.humancellatlas.org/).
All data from the HCA that we are using can be found in the S3 bucket, `s3://ccdl-scpca-data/human_cell_atlas_data`. 
The following files all contain information related to the test datasets. 

1. `hca-project-metadata.tsv`: This file contains information about each of the projects that are being used for testing from the HCA. 
Each row in this file corresponds to a project and contains the following columns: 

| column_id         | contents                                                           |
|-------------------|--------------------------------------------------------------------|
| `project_name`    | The shorthand project name assigned by the HCA                     |
| `source_id`       | Unique ID associated with the project                              |
| `tissue_group`    | Tissue group the project belongs to (e.g. blood, brain, or kidney) |
| `files_directory` | files directory on S3                                              |
| `metadata_filename`| name of metadata file obtained from the HCA                       |

2. `hca-library-metadata.tsv` This file contains information about each library from every project that is being used as a test dataset from the HCA. 
Each row in this file corresponds to a library and contains the following columns: 

| column_id         | contents                                                           |
|-------------------|--------------------------------------------------------------------|
| `sample_biomaterial_id`    | Unique ID associated with the individual sample           |
| `library_biomaterial_id`   | Unique ID associated with the individual library that was sequenced |
| `bundle_uuid`    | UUID for the individual folder containing each loom file |
| `project_name` | The shorthand project name assigned by the HCA                        |
| `source_id`| Unique ID associated with the project                                     |
| `tissue_group`    | Tissue group the project belongs to (e.g. blood, brain, or kidney) |
| `technology`       |  Sequencing/library technology used (e.g. 10Xv2, 10Xv3, etc.)          |
| `seq_unit`    | Sequencing unit (cell or nucleus)                         |
| `diagnosis`       | Indicates if the sample came from diseasead or normal tissue       |
| `organ`    | Specified tissue by the HCA where the sample was obtained from            |
| `organ_part` | Specified tissue region by the HCA where the sample was obtained from                                              |
| `selected_cell_types`| Identifies the group of cells selected for prior to sequencing, otherwise NA |
| `s3_files_dir`    | files directory on S3                                              |    
| `loom_file`       | loom file name in the format `tissue_group/project_name/bundle_uuid/filename` |

3. `hca-processed-libraries.tsv`: This file contains the list of libraries from each project that are being used for testing data integration. 
This file is used as input to the script, `scripts/00-convert-loom.R`, used for converting loom to `SingleCellExperiment` objects and saving those objects to S3. 
