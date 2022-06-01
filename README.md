# sc-data-integration

This repo contains files, scripts, and analysis related to exploring integration of single-cell and single-nuclei data. 

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
| `library_biomaterial_id`   | Unique ID associated with the individual library that was sequenced                             |
| `bundle_uuid`    | UUID for the individual folder containing each loom file |
| `project_name` | The shorthand project name assigned by the HCA                        |
| `source_id`| Unique ID associated with the project                                     |
| `tissue_group`    | Tissue group the project belongs to (e.g. blood, brain, or kidney) |
| `technology`       | Unique ID associated with the project                             |
| `seq_unit`    | The shorthand project name assigned by the HCA                         |
| `diagnosis`       | Unique ID associated with the project                              |
| `organ`    | Specified tissue by the HCA where the sample was obtained from            |
| `organ_part` | Specified tissue region by the HCA where the sample was obtained from                                              |
| `selected_cell_types`| Identifies the group of cells selected for prior to sequencing, otherwise NA                       |
| `s3_files_dir`    | files directory on S3                                              |    
| `loom_file`       | loom file name in the format `tissue_group/project_name/bundle_uuid/filename`                            |
