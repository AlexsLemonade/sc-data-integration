#!/usr/bin/env python3
"""
Script to test integration with scanorama
"""

import os
import git

project_root = git.Repo('.', search_parent_directories=True).working_dir

import anndata as adata
from integrate_scanorama import integrate_scanorama

anndata_dir = os.path.join(project_root, 
                           "results",
                           "human_cell_atlas", 
                           "anndata")

# build merged anndata file paths 
merged_anndata_dir = os.path.join(anndata_dir,
                                  "merged_anndata_objects")
merged_adata_file = os.path.join(merged_anndata_dir, "1M_Immune_Cells_anndata.h5")

# build integrated anndata file paths for output
integrated_adata_dir = os.path.join(anndata_dir,
                                    "integrated_scanorama_objects")
scanorama_integrated_adata_file = os.path.join(integrated_adata_dir,
                                    "1M_Immune_Cells_scanorama_integrated.h5")

# check if output directory exists 
if os.path.isdir(integrated_adata_dir)==False:
    os.makedirs(integrated_adata_dir)


# read in merged anndata object
merged_adata = adata.read_h5ad(merged_adata_file)

# integrate anndata with scanorama
scanorama_integrated_adata = integrate_scanorama(merged_adata, seed = 2022)

# write anndata to h5
scanorama_integrated_adata.write(filename = scanorama_integrated_adata_file)