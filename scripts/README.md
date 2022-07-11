# Scripts for single-cell data integration

This folder holds the scripts that have been used for data integration and processing of single-cell libraries and datasets.

## Obtaining SingleCellExperiment objects 

The `00-obtain-sce.R` script is specifically for working with the test datasets from the Human Cell Atlas Data portal. 
Data downloaded directly from the portal are available as `SingleCellLoomExperiment` objects as a loom file, but for all downstream processing, `SingleCellExperiment` objects as RDS files will be required as input. 
If the `SingleCellExperiment` objects have already been created by another Data Lab staff member, this script can be used to sync `SingleCellExperiment` objects from the S3 bucket to a local environment.
If working with libraries for the first time that don't yet have `SingleCellExperiment` objects that exist on S3, they will be created by converting the `loom` files and then syncing the `SingleCellExperiment` object to S3. 

All libraries that are being used for testing should be listed in `sample-info/hca-processed-libraries.tsv`. 
Before starting analysis, to copy the `SingleCellExperiment` objects for libraries listed in `scripts/hca-processed-libraries.tsv` locally, use the `--copy_s3` flag. 
This will copy the `SingleCellExperiment` objects to the `data/human_cell_atlas/sce` folder where libraries will be nested by tissue group and project name. 

```R
Rscript 00-obtain-sce.R --copy_s3
```

To convert `loom` files to `SingleCellExperiment` objects as RDS files for the first time run the script with the default settings. 

```R
Rscript 00-obtain-sce.R
```

If changes have been made to the way files are converted and you would like to re-convert the `loom` files, use the `--overwrite` option. 

```R
Rscript --obtain-sce.R --overwrite
```

## Running HCA test datasets through scpca-downstream-analyses 

Before libraries can be integrated, data downloaded from the Human Cell Atlas Data portal will need to be filtered to remove empty droplets and low quality cells and undergo normalization. 
The `01-run-downstream-analyses.sh` script is a bash script that should be used to both remove empty drops and run the [core workflow as part of `scpca-downstream-analyses`](https://github.com/AlexsLemonade/scpca-downstream-analyses#the-core-downstream-analyses-workflow) for each library listed in `sample-info/hca-processed-libraries.tsv`. 

The bash script first runs the Rscript `utils/preprocess-sce.R`, which will create a filtered `SingleCellExperiment` object for each library. 
The filtered results can be found in `results/human_cell_atlas/filtered_sce`. 
This preprocessing script also creates a metadata file, `sample-info/hca-downstream-metadata.tsv`, that is required to run `scpca-downstream-analyses` containing the following columns: 

- `sample ID`: Equivalent to the `sample_biomaterial_id`
- `library ID`: Equivalent to the `library_biomaterial_id`
- `filtering_method`: What type of filtering to use to remove low quality cells. 
This is set to use [`miQC`](https://rdrr.io/github/greenelab/miQC/man/filterCells.html), but if `miQC` fails, a set of manual filtering thresholds will be used instead.
- `filepath`: The full path to the filtered `SingleCellExperiment` object as an RDS file.
This path is specific to your local machine. 

Following creation of the metadata, `scpca-downstream-analyses` is run for each library. 
This produces a filtered, normalized `SingleCellExperiment` object and a summary html report showing the filtering and clustering results for that library. 
Results will be stored in `results/human_cell_atlas/scpca-downstream-analyses`.

Running the script with default settings will produce both the empty drops filtered output and the output from `scpca-downstream-analyses`.
You will need to first clone the [`scpca-downstream-repo`](https://github.com/AlexsLemonade/scpca-downstream-analyses) and provide the full path to the location of the repo on your local computer.

```R
bash 01-run-downstream-analyses.sh --downstream_repo <path to location of scpca-downstream-repo>
```

If desired, an S3 bucket can be provided to ensure that the results from `scpca-downstream-analyses` are copied to S3 for other Data Lab staff members to use. 

```R 
bash 01-run-downstream-analyses.sh \
  --downstream_repo <path to location of scpca-downstream-repo> \
  --s3_bucket "s3://sc-data-integration/human_cell_atlas_results"
```
