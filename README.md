# sc-data-integration
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Managing R packages with `renv`](#managing-r-packages-with-renv)
- [Metadata](#metadata)
- [Shared data files](#shared-data-files)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

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
This file is used as input to the script, `scripts/00-obtain-sce.R`, used for converting loom to `SingleCellExperiment` objects and saving those objects to S3. 

## Shared data files 

Data for this project is stored in the S3 bucket, `s3://sc-data-integration`. 
The following data can be found in the above S3 bucket within the `human_cell_atlas_data` folder: 

- The `loom` folder contains the original loom files downloaded from the Human Cell Atlas data portal for each test dataset. 
Here loom files are nested by `tissue_group`, `project_name`, and `bundle_uuid`. 
- The `sce` folder contains the unfiltered `SingleCellExperiment` objects saved as RDS files.
These `SingleCellExperiment` objects have been converted from the loom files using the `00-obtain-sce.R` script in the `scripts` directory in this repo.
Here RDS files are nested by `tissue_group` and `project_name`.

A separate `reference-files` folder contains any reference files needed for processing dataset, such as the gtf file needed to generate the mitochondrial gene list found in the `reference-files` folder in the repository. 

In order to access these files, you must have credentials setup for AWS. 
Additionally, some of the scripts in this repository require use of AWS command line tools. 
We have previously written up detailed instructions on [installing the AWS command line tools and configuring your credentials](https://github.com/AlexsLemonade/alsf-scpca#aws) that can be used as a reference.

After AWS command line tools have been set up, the `SingleCellExperiment` objects found in `s3://sc-data-integration/human_cell_atlas_data/sce` can be copied to your local computer by running the `00-obtain-sce.R` script with the `--copy_s3` flag.

```
Rscript scripts/00-obtain-sce.R --copy_s3
```

This will copy any `SingleCellExperiment` objects for libraries listed in `hca-processed-libraries.tsv` that have already been converted from loom files.
If any libraries listed in `hca-processed-libraries.tsv` do not have corresponding `SingleCellExperiment` objects, running the `00-obtain-sce.R` will also convert those loom files.
