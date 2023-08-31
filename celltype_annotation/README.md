# Cell type annotation

This folder holds analysis and scripts developed for exploring cell type annotation methods.

## Analysis

Inside the analysis folder are a set of template notebooks that have been used for exploring cell type annotations using different methods.

1. `SingleR-non-immune-ref-comparison.Rmd`: This notebook was used to compare cell type annotations obtained from `SingleR`.
We compared the annotations when using references with or without immune cell populations.
In particular, this notebook is specifically useful for looking at cell type annotation in blood tumors.

2. `SingleR-combining-refs.Rmd`: This notebook specifically compares using multiple `celldex` references to using a single `celldex` reference to annotate libraries with `SingleR`.
Again, this notebook specifically is for blood tumors.

3. `marker-gene-celltype-exploration.Rmd`: This notebook looks at using other marker gene based methods, specifically `scType` and `CellAssign`.

4. `01-cell-assign-sarcoma.Rmd`: Here we specifically look at using `CellAssign` with Rhabdomyosarcoma samples from `SCPCP000005`.
This analysis compares the assignments using `CellAssign` to the submitter annotations.

5. `02-cell-assign-sarcoma-marker-genes.Rmd`: This notebook also compare annotations using `CellAssign` to manual annotations cell type annotations in Rhabdomyosarcoma samples from `SCPCP000005`.
However instead of using public references for marker gene sets, we use `scran::FindMarkers()` to build a custom marker gene table for the specified library.

6. `singler-rms-comparison.Rmd`: This notebook compares annotations from `SingleR` to manual annotations in Rhabdomyosarcoma samples from `SCPCP000005`.

7. `03-cell-assign-combined-markers.Rmd`: This notebook is follow-up to the exploration in `01-cell-assign-sarcoma.Rmd`.
In the notebook, cell type annotations using different sets of marker genes from `PanglaoDB` are compared to submitter annotations.

8. `copykat-rms-exploration.Rmd`: This notebook explores the usage of `CopyKAT` for identifying tumor vs. normal cells in Rhabdomyosarcoma samples from `SCPCP000005`.

9. `04-cell-assign-delta-median.Rmd`: This notebook examines ways to evaluate confidence in the cell type assignments that are obtained from using `CellAssign`. 
The notebook uses the cell types annotated in `marker-gene-celltype-exploration.Rmd` ,`02-cell-assign-sarcoma-marker-genes.Rmd`, and `03-cell-assign-combined-markers.Rmd`.

## Scripts

This folder contains any scripts used for cell type annotation.

1. `cell-assign.py`: This script can be used to perform cell type annotation using `CellAssign`.
The output will be a table of predictions.
2. `remove-cell-refs.py`: This script is used to remove a set of cell types from a `celldex` reference.
This requires a txt file with a list of cell types to remove.
The default is the `immune_celltypes.txt` file in the root directory of this repo.
3. `render-SingleR-reports.R`: This script is used to render the `SingleR-non-immune-ref-comparison.Rmd` report across an entire project.

## Utils

This folder holds helper functions that were used across multiple analysis templates.

1. `cellassign-helper-functions.R`: This script contains any helper functions used for analyzing cell type annotations obtained from `CellAssign`.
2. `singler-helper-functions.R`: his script contains any helper functions used for analyzing cell type annotations obtained from `SingleR`.
