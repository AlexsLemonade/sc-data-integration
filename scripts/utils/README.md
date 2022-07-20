# Utilities for single-cell data integration

This directory holds the scripts with functions that are used for data integration and processing of single-cell libraries and datasets.


- The script `preprocess-sce.R` is used to filter SCE objects using [`scpcaTools::filter_counts`](https://github.com/AlexsLemonade/scpcaTools/blob/main/R/filter_counts.R) and generates a metadata file for use by the [`scpca-downstream-analyses`](https://github.com/AlexsLemonade/scpca-downstream-analyses/) pipeline. 
To run this script followed by processing SCE through `scpca-downstream-analyses` to complete filtering, normalization, and clustering use the `01-run-downstream-analyses.sh` script in the `scripts` directory. 
- The script `convert-sce-to-anndata.R` converts SCE objects stored as RDS files to [`anndata`](https://anndata.readthedocs.io/en/latest/) objects stored as HDF5 files for input to the integration method [`Scanorama`](https://github.com/brianhie/scanorama) within the Python package [`scanpy`](https://github.com/scverse/scanpy).


- There are several scripts that contain functions to perform integration with various methods:
    - `integrate-harmony.R` contains the function to perform single-cell library integration using the [`harmony`](https://github.com/immunogenomics/harmony) R package function `harmonyMatrix()`.
    - `integrate-fastMNN.R` contains the function to perform single-cell library integration using the [`batchelor`](https://bioconductor.org/packages/devel/bioc/html/batchelor.html) R package function `fastMNN()`.

- The script `integration-helpers.R` contains several functions that support integration:
    - `perform_hvg_selection()` identifies high-variance genes to use during integration.
    - `perform_dim_reduction()` calculates PCA (if specified) and UMAP embeddings from a combined SCE object.
    - `combine_sce_objects()` merges one or more SCE objects that should be integrated.

- The script `test_integration-functions.R` is used internally to lightly test functionality of new features. 
Note that this script does _not_ facilitate true unit testing/regression testing. 

