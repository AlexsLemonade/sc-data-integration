# sc-data-integration
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Environment setup](#environment-setup)
  - [Managing R packages with `renv`](#managing-r-packages-with-renv)
  - [Snakemake/conda setup](#snakemakeconda-setup)
- [Data used for benchmarking integration](#data-used-for-benchmarking-integration)
  - [Metadata](#metadata)
  - [Data files present on S3](#data-files-present-on-s3)
  - [Processed SingleCellExperiment objects to use for data integration](#processed-singlecellexperiment-objects-to-use-for-data-integration)
- [Running the integration workflow](#running-the-integration-workflow)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

This repo contains files, scripts, and analysis related to exploring integration of single-cell and single-nuclei data.

## Environment setup

### Managing R packages with `renv`

Package dependencies for the analysis workflows in this repository are managed using [`renv`](https://rstudio.github.io/renv/index.html).
For `renv` to work as intended, you'll need to work within the `sc-data-integration.Rproj` project in RStudio.
You may need to run `renv::restore()` upon opening the project to ensure the `renv.lock` file is synced with the project library.

Each time you install or use new packages, you will want to run `renv::snapshot()` to update the `renv.lock` file with any added package and dependencies necessary to run the analyses and scripts in this repo.

If there are dependencies you want to include that are not captured automatically by `renv::snapshot()`, add them to `components/dependencies.R` with a call to `library()` and an explanatory comment.
For example, if `dplyr` were recommended but not required by a package and you wanted to make sure to include it in the lockfile, you would add `library(dplyr)` to `components/dependencies.R`.
Then rerun `renv::snapshot()`.

### Snakemake/conda setup

The main workflow for the integration scripts is written with Snakemake, which will handle most dependencies internally, including the `renv` environment.

You will need the latest version of `snakemake` and the `peppy` python package.
The easiest way to install these is with `conda` and/or `mamba`, which you will want to set up to use the `bioconda` and `conda-forge` channels using the following series of commands:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

You can then install `snakemake` and `peppy` into your preferred environment with:

```
mamba install snakemake peppy
```

(Use `conda install` if you do not have `mamba` installed.)



Python-based environments will be built automatically by Snakemake when the workflow is run, but the environment for R should be built before running the workflow.
To create or update the necessary environment for the R scripts, which includes an isolated version of `R`, `pandoc`, and the `renv` package installation, run the following command from the base of the repository:

```
bash setup_envs.sh
```
This script will use Snakemake to install all necessary components for the workflow in an isoloated Conda enviroment.
If you are on an Apple Silicon (M1/M2/Arm) Mac, this should properly handle setting up R to use an Intel-based build for compatibility with Bioconductor packages.

This installation may take up to an hour, as all of the R packages will likely have to be compiled from scratch.
However, this should be a one-time cost, and ensures that you have all of the tools for the workflow installed and ready.

To use the environment you have just created, you will need to run Snakemake with the `--use-conda` flag each time.

If there are updates to the `renv.lock` file only, those can be applied with the following command (on any system):

```
snakemake --use-conda -c2 build_renv
```

## Data used for benchmarking integration

For exploring data integration, we used test datasets that were obtained from the [Human Cell Atlas (HCA) Data Portal](https://data.humancellatlas.org/), the [Single-cell Pediatric Cancer Atlas(ScPCA) Portal](scpca.alexslemonade.org), and simulated single-cell data published in [Luecken _et al.,_ (2022)](https://doi.org/10.1038/s41592-021-01336-8).

All data from the HCA that we are using can be found in the private S3 bucket, `s3://sc-data-integration/human_cell_atlas_data`.
The simulated data can be downloaded directly from [figshare](https://doi.org/10.6084/m9.figshare.12420968).

All gene expression data used for benchmarking is stored in the private S3 bucket, `s3://sc-data-integration`.
In order to access these files, you must be a Data Lab staff member and have credentials setup for AWS.

### Metadata

Inside the `sample-info` folder is metadata related to datasets used for testing data integration.

1. `<project_name>-project-metadata.tsv`: This file contains information about each of the projects that are being used for testing integration from a given area (e.g., HCA, simulated, ScPCA).
Each row in this file corresponds to a project, dataset, or group of libraries that should be integrated together.
All `project-metadata.tsv` files must contain a `project_name` column, but may also contain other relevant project information such as the following:

| column_id         | contents                                                           |
|-------------------|--------------------------------------------------------------------|
| `project_name`    | The shorthand project name                     |
| `source_id`       | Unique ID associated with the project                              |
| `tissue_group`    | Tissue group the project belongs to (e.g. blood, brain, or kidney) |
| `files_directory` | files directory on S3                                              |
| `metadata_filename`| name of metadata file obtained from the HCA                       |
| `celltype_filename` | file name corresponding to file containing cell type information as found on HCA                                              |
| `celltype_filetype`| format of cell type file availble on HCA                       |

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

3. `<project_name>-processed-libraries.tsv`: This file contains the list of libraries from each project that are being used for testing data integration.
This metadata file is required for most scripts, including running `scpca-downstream-analyses` using `01-run-downstream-analyses.sh` and for running the integration workflow.
This file must contain the following columns, but may also contain additional columns related to a given dataset:

| column_id         | contents                                                           |
|-------------------|--------------------------------------------------------------------|
| `sample_biomaterial_id`    | Unique ID associated with the individual sample           |
| `library_biomaterial_id`   | Unique ID associated with the individual library that was sequenced |
| `project_name`    | The shorthand project name                 |
| `integration_input_dir` | The directory containing the `SingleCellExperiment` objects to be used as input to the data integration snakemake workflow |


4. `hca-celltype-info.tsv`: *This file is not available on the repo and is stored in the private S3 bucket, `s3://sc-data-integration/sample-info`.*
This file contains all available cell type information for projects listed in `hca-project-metadata.tsv`.
This file was created using the `scripts/00a-reformat-celltype-info.R` which takes as input the cell type information available for each project from the Human Cell Atlas Data Portal.
The cell type information for each project, in its original format, can be stored in `s3://sc-data-integration/human_cell_atlas_data/celltype`.
Each row corresponds to a single cell and contain the following information:

| column_id         | contents                                                           |
|-------------------|--------------------------------------------------------------------|
| `sample_biomaterial_id`    | Unique ID associated with the individual sample           |
| `library_biomaterial_id`   | Unique ID associated with the individual library that was sequenced |
| `project` | The shorthand project name assigned by the HCA                        |
| `barcode` | The unique cell barcode                       |
| `celltype` | The assigned cell type for a given barcode, obtained from cell type data stored in `s3://sc-data-integration/human_cell_atlas_data/celltype`                        |


### Data files present on S3

All data and intermediate files are stored in the private S3 bucket, `s3://sc-data-integration`.
The following data can be found in the above S3 bucket within the `human_cell_atlas_data` folder:

- The `loom` folder contains the original loom files downloaded from the Human Cell Atlas data portal for each test dataset.
Here loom files are nested by `tissue_group`, `project_name`, and `bundle_uuid`.
- The `sce` folder contains the unfiltered `SingleCellExperiment` objects saved as RDS files.
These `SingleCellExperiment` objects have been converted from the loom files using the `00-obtain-sce.R` script in the `scripts` directory in this repo.
Here RDS files are nested by `tissue_group` and `project_name`.

The following data can be found in the S3 bucket within the `scib_simulated_data` folder:

- The `hdf5` folder contains the original `hdf5` files for simulated data obtained from [figshare](https://doi.org/10.6084/m9.figshare.12420968).
- The `sce` folder contains the individual `SingleCellExperiment` objects stored as `rds` files after running `scripts/00b-obtain-sim-sce.R` and `scripts/00c-create-sim1-subsets.R`

A separate `reference-files` folder contains any reference files needed for processing dataset, such as the gtf file needed to generate the mitochondrial gene list found in the `reference-files` folder in the repository.

In order to access these files, you must be a Data Lab staff member and have credentials setup for AWS.
Additionally, some of the scripts in this repository require use of AWS command line tools.
We have previously written up detailed instructions on [installing the AWS command line tools and configuring your credentials](https://github.com/AlexsLemonade/alsf-scpca#aws) that can be used as a reference.

After AWS command line tools have been set up, the `SingleCellExperiment` objects found in `s3://sc-data-integration/human_cell_atlas_data/sce` can be copied to your local computer by running the `00-obtain-sce.R` script with the `--copy_s3` flag.

```
Rscript scripts/00-obtain-sce.R --copy_s3
```

This will copy any `SingleCellExperiment` objects for libraries listed in `hca-processed-libraries.tsv` that have already been converted from loom files.
If any libraries listed in `hca-processed-libraries.tsv` do not have corresponding `SingleCellExperiment` objects, running the `00-obtain-sce.R` will also convert those loom files.

### Processed SingleCellExperiment objects to use for data integration

The `human_cell_atlas_results/scpca-downstream-analyses` folder contains all processed `SingleCellExperiment` objects and the output from [running the core workflow in `scpca-downstream-analyses`](https://github.com/AlexsLemonade/scpca-downstream-analyses/blob/main/Snakefile).
Within this folder each library that has been processed has its own folder that contains both the processed `SingleCellExperiment` object and an html summary report.
The `SingleCellExperiment` objects in this folder have both been filtered to remove empty droplets and run through `scpca-downstream-analyses` using the `scripts/01-run-downstream-analyses.sh` script.
This means they contain a `logcounts` assay with the normalized counts matrix, both PCA and UMAP embeddings, and clustering assignments that can be found in the `louvain_10` column of the `colData`.
The `SingleCellExperiment` objects present in `human_cell_atlas_results/scpca-downstream-analyses` should be the objects used as input for integration methods.

These files were produced and synced to S3 using the following script:

**Note:** To run the below script, you must have available in your path R (v4.1.2), [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba) and [`pandoc`](https://pandoc.org/installing.html#macos).
`pandoc` must be version 1.12.3 or higher, which can be checked using the `pandoc -v` command.

```
bash scripts/01-run-downstream-analyses.sh \
  --downstream_repo <full path to scpca-downstream-analyses-repo> \
  --s3_bucket "s3://sc-data-integration/human_cell_atlas_results/scpca-downstream-analyses"
```

**Note:** If you wish to run ScPCA data through the integration workflow, rather than HCA data, please see the [special guidelines for preparing ScPCA data](./scripts/README.md#preparing-scpca-test-datasets-for-integration).

## Running the integration workflow

To run the integration workflow, invoke `snakemake` from the `sc-data-integration` directory with the following command:

```
snakemake -c4 --use-conda
```

You can adjust the number of cores used by adjusting the `-c4` flag with however many cores you want to use where the given number represents the number of desired cores (here, 4).
Note that you will want to have [set up the R conda environment already](#snakemakeconda-setup), especially if you are on an Apple Silicon Mac.

To run the workflow for development, you may wish to specify the `config-test.yaml` file, which will only run one project through the pipeline to save time:

```
snakemake -c4 --use-conda --configfile config-test.yaml
```

Finally, to run the `scib_simulated` data through the pipeline, use:

```
snakemake -c4 --use-conda --configfile config-scib_simulated.yaml
```

