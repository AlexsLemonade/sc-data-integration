# Scripts for single-cell data integration

This directory holds the scripts that have been used for data integration and processing of single-cell libraries and datasets.
The `utils/` directory contains functions used within these scripts and is separately documented in [`utils/README.md`](./utils/README.md).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Obtaining `SingleCellExperiment` objects](#obtaining-singlecellexperiment-objects)
- [Re-formatting cell type information](#re-formatting-cell-type-information)
  - [Special considerations for cell type information](#special-considerations-for-cell-type-information)
- [Running HCA test datasets through `scpca-downstream-analyses`](#running-hca-test-datasets-through-scpca-downstream-analyses)
  - [Generating the mitochondrial gene list](#generating-the-mitochondrial-gene-list)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Obtaining `SingleCellExperiment` objects

The `00-obtain-sce.R` script is specifically for working with the test datasets from the [Human Cell Atlas Data portal](https://data.humancellatlas.org/).
Data downloaded directly from the portal are available as `SingleCellLoomExperiment` objects as `loom` files, but for all downstream processing, `SingleCellExperiment` objects as RDS files will be required as input.
If the `SingleCellExperiment` objects have already been created by another Data Lab staff member, this script can be used to sync `SingleCellExperiment` objects from the S3 bucket to a local environment.
If working with libraries for the first time that don't yet have `SingleCellExperiment` objects that exist on S3, they will be created by converting the `loom` files and then syncing the `SingleCellExperiment` objects to S3.

All libraries that are being used for testing should be listed in [`sample-info/hca-processed-libraries.tsv`](../sample-info/hca-processed-libraries.tsv).
Before starting analysis, copy the `SingleCellExperiment` objects for libraries listed in `scripts/hca-processed-libraries.tsv` locally using the `--copy_s3` flag.
This will copy the `SingleCellExperiment` objects to the `data/human_cell_atlas/sce` directory where libraries will be nested by tissue group and project name.
This will also copy all of the available cell type information in its original format from the Human Cell Atlas Data Portal, and save it to `data/human_cell_atlas/celltype`, where files will be nested by tissue group and project name.

```R
Rscript 00-obtain-sce.R --copy_s3
```

To convert `loom` files to `SingleCellExperiment` objects as RDS files for the first time, run the script with the default settings.
This will only perform the conversion for new libraries added to `sample-info/hca-processed-libraries.tsv` that have not been converted previously and are not present on S3.

```R
Rscript 00-obtain-sce.R
```

If changes have been made to the way files are converted and you would like to re-convert the `loom` files, use the `--overwrite` flag.

```R
Rscript 00-obtain-sce.R --overwrite
```

## Re-formatting cell type information

The cell type information for each project obtained from the Human Cell Atlas Data portal is available in its own unique format.
No two cell type files are alike!

The original, unedited files can be found at `s3://sc-data-integration/human_cell_atlas_data/celltype`, and copied locally using the `00-obtain-sce.R` script with the `--copy_s3` flag as described above.
After these files are available locally in the `data/human_cell_atlas/celltype` folder, they can be reformatted into a single file containing the cell types for all cells across all libraries and projects ([`sample-info/hca-celltype-info.tsv`](../sample-info/hca-celltype-info.tsv)).
The `00a-reformat-celltype-info.R` script reformats the files from all projects with available cell type information and then combines into a single individual file.

### Special considerations for cell type information

A few of the projects required some extra changes that were performed prior to running the `00a-reformat-celltype-info.R` script.

1. `HumanBrainSubstantiaNigra` - For this project, celltype information was available in an excel spreadsheet containing two captioned tables.
The table containing the cortex samples was copied and saved as a tsv file, `cortex_only_41467_2020_17876_MOESM5_ESM.tsv` which was read in directly to the reformatting script.
Additionally, the `library_biomaterial_id`corresponded to an SRA identifier, while the sample names found in the provided excel sheet were not SRA identifiers.
The SRA identifier was searched and matched to the corresponding sample identifier used by the submitter.
2. `Oligodendrocyte_MS` - The cell type information was available as a txt file which was first converted to a tsv file.
Additionally, the `library_biomaterial_id`corresponded to an SRA identifier, while the sample names found in the provided txt file were not SRA identifiers.
The SRA identifier was searched and matched to the corresponding sample identifier used by the submitter.
3. `MultipleSclerosisLineageDiversity`- For this project, celltype information was found in a multi-sheet excel file where each file contained a different supplemental table.
The celltype information was only found in supplemental table 2 which was copied and saved as a tsv file, `cell_metadata_EMS83370-supplement-2.tsv` which was read in directly to the reformatting script.
Additionally, the `library_biomaterial_id`corresponded to an SRA identifier, while the sample names found in the provided excel sheet were not SRA identifiers.
The SRA identifier was searched and matched to the corresponding sample identifier used by the submitter.
4. `FetalLiverHaematopoiesis` - The library ID's provided in the cell type file were slightly different then the `library_biomaterial_id` available on the HCA portal.
The ID submitted by the user in the cell type file was mapped to the `library_biomaterial_id` available in the HCA portal.


## Running HCA test datasets through `scpca-downstream-analyses`

Before libraries can be integrated, data downloaded from the Human Cell Atlas Data portal will need to be filtered to remove empty droplets and low quality cells and undergo normalization.
The `01-run-downstream-analyses.sh` script is a bash script that should be used to both remove empty drops and run the [core workflow from `scpca-downstream-analyses`](https://github.com/AlexsLemonade/scpca-downstream-analyses#the-core-downstream-analyses-workflow) for each library listed in `sample-info/hca-processed-libraries.tsv`.

The bash script first runs the R script `utils/preprocess-sce.R`, which will create a filtered `SingleCellExperiment` object for each library.
The filtered results can be found in `results/human_cell_atlas/filtered_sce`.
This preprocessing script also creates a metadata file, `sample-info/hca-downstream-metadata.tsv`, that is required to run `scpca-downstream-analyses`, and contains the following columns:

- `sample_id`: Equivalent to the `sample_biomaterial_id`
- `library_id`: Equivalent to the `library_biomaterial_id`
- `filtering_method`: What type of filtering to use to remove low quality cells.
This is set to use [`miQC`](https://rdrr.io/github/greenelab/miQC/man/filterCells.html), but if `miQC` fails, a set of manual filtering thresholds will be used instead.
- `filepath`: The full path to the filtered `SingleCellExperiment` object as an RDS file.
This path is specific to your local machine.

Following creation of the metadata, the `scpca-downstream-analyses` workflow is run for each library.
This produces a filtered, normalized `SingleCellExperiment` object and a summary html report showing the filtering and clustering results for that library.
Results will be stored in `results/human_cell_atlas/scpca-downstream-analyses`.

Running the script with default settings will produce both the empty drops filtered output and the output from `scpca-downstream-analyses`.

To run `01-run-downstream-analyses.sh`, you will need to first clone the [`scpca-downstream-repo`](https://github.com/AlexsLemonade/scpca-downstream-analyses) and provide the full path to the location of the repo on your local computer.
This script requires a mitochondrial gene list to use for calculating the mitochondrial reads present in each cell.
For datasets downloaded from the Human Cell Atlas Data portal, the mitochondrial gene list has already been created using `generate-mito-reference.R` (see [instructions below](#generating-the-mitochondrial-gene-list)) and can be found in [`reference-files/gencode.v27.mitogenes.txt`](../reference-files/gencode.v27.mitogenes.txt).
This file is the default file that is used in running `01-run-downstream-analyses.sh`.
**Note:** You must have available in your path R (v4.1.2),  [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba) and [`pandoc`](https://pandoc.org/installing.html#macos).
`pandoc` must be version 1.12.3 or higher, which can be checked using the `pandoc -v` command.

```sh
bash 01-run-downstream-analyses.sh --downstream_repo <path to location of scpca-downstream-repo>
```

If desired, an S3 bucket link can be provided to ensure that the results from `scpca-downstream-analyses` are copied to S3 for other Data Lab staff members to use.

```sh
bash 01-run-downstream-analyses.sh \
  --downstream_repo <path to location of scpca-downstream-repo> \
  --s3_bucket "s3://sc-data-integration/human_cell_atlas_results"
```

By default, filtering of empty droplets will not be repeated if filtered `SingleCellExperiment` objects are already present locally.
To overwrite existing filtered files, use `--repeat_filtering yes` at the command line.

```sh
bash 01-run-downstream-analyses.sh \
  --downstream_repo <path to location of scpca-downstream-repo> \
  --repeat_filtering yes
```

### Generating the mitochondrial gene list

To process libraries through the core workflow in `scpca-downstream-analyses`, a text file containing a list of mitochondrial genes found in the reference transcriptome used for alignment is required as input.
Within this file, each row must contain a unique ensembl gene identifier corresponding to a mitochondrial gene.

To generate this list for the reference used, use the `generate-mito-reference.R` script.
Running the script with the default configuration will re-create the `reference-files/gencode.v27.mitogenes.txt` mitochondrial gene list.

```R
Rscript generate-mito-reference.R
```

To generate a mitochondrial gene list for a different reference, you will need to obtain the unzipped gtf file corresponding to the reference genome or transcriptome of interest.
The gtf file should be stored in the `reference-files` directory before running the script.

```R
Rscript generate-mito-reference.R \
  --gtf_file <filename for gtf file> \
  --mito_output <full path to save mitochondrial gene list>
```
