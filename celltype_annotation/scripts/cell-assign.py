#!/usr/bin/env python3

## basic script to use for testing cellassign with
# followed tutorial from scvi cellassign
# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/cellassign_tutorial.html
# used Follicular lympyhoma reference on zenodo
# https://zenodo.org/record/3372746#.ZDcGQ-zMJA0


import os
import anndata as adata
import scvi
from scvi.external import CellAssign
import pandas as pd
import numpy as np

# read in public references as marker gene tables
fl_file = os.path.join("..", "references", "FL_celltype_ensembl.csv")
fl_mtx = pd.read_csv(fl_file, index_col=0)

fl_noBcells_file = os.path.join("..", "references", "FL_celltype_ensembl-noBcells.csv")
fl_noBcells_mtx = pd.read_csv(fl_noBcells_file, index_col=0)

# file path to annotated sce
anndata_dir = os.path.join("..", "data", "anndata")
annotated_adata_file = os.path.join(anndata_dir, "SCPCL000295_anndata.hdf5")
annotated_adata = adata.read_h5ad(annotated_adata_file)

# index adata but use row data to do it
subset_adata = annotated_adata[:, fl_mtx.index].copy()
subset_adata.X = subset_adata.X.tocsr()

# renormalize/ calcualate size factor on subset adata
lib_size = annotated_adata.X.sum(1)
subset_adata.obs["size_factor"] = lib_size / np.mean(lib_size)

# small function for get cell type predictions given a marker genes matrix
def get_predictions(anndata_object, marker_genes_mtx):
    scvi.external.CellAssign.setup_anndata(anndata_object, size_factor_key="size_factor")
    model = CellAssign(anndata_object, marker_genes_mtx)
    model.train()
    predictions = model.predict()
    return(predictions)

# model cell types with all cells
all_cells_predictions = get_predictions(anndata_object = subset_adata,
                                        marker_genes_mtx = fl_mtx)


# model cell types with no B cells
no_bcells_predictions = get_predictions(anndata_object = subset_adata,
                                        marker_genes_mtx = fl_noBcells_mtx)

# add predictions to annotated adata
annotated_adata.obs["fl_all_prediction"] = all_cells_predictions.idxmax(axis=1).values
annotated_adata.obs["fl_noBcells_prediction"] = no_bcells_predictions.idxmax(axis=1).values

# write out modified adata
output_anndata = os.path.join(anndata_dir, "SCPCL000295_cellassign.hdf5")
annotated_adata.write(filename = output_anndata)

# write out predictions as tsv
all_cells_predictions.to_csv(os.path.join(anndata_dir, "all_cells_predictions.tsv"), sep = "\t", index = False)
no_bcells_predictions.to_csv(os.path.join(anndata_dir, "no_bcells_predictions.tsv"), sep = "\t", index = False)

