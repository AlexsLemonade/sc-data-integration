#!/usr/bin/env python3
"""
Script to perform integration on a merged AnnData object using Scanorama

Takes
"""

import os
project_root = git.Repo('.', search_parent_directories=True).working_dir
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


import git
import anndata as adata
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--input_anndata',
                    default = merged_adata_file,
                    help = 'Path to HDF5 file with merged AnnData object to integrate')
parser.add_argument('-i', '--output_anndata',
                    default = scanorama_integrated_adata_file,
                    help = 'Path to HDF5 file to save the integrated AnnData object')
parser.add_argument('-b', '--batch_column',
                    default = 'batch',
                    help = 'The name of the column in `anndata.obs` that indicates the batches for each cell, " \
                            " typically this corresponds to the library id.')
parser.add_argument('-s', '--seed',
                    default = None,
                    help = 'Random seed to set prior to scanorama.')

from utils.integrate_scanorama import integrate_scanorama

# check that input file has correct extension

# check that batch column in anndata? or since function has check we don't need it?

# read in merged anndata object
merged_adata = adata.read_h5ad(args.input_anndata)

# integrate anndata with scanorama
scanorama_integrated_adata = integrate_scanorama(merged_adata,
                                                 batch_column = args.batch_column,
                                                 seed = args.seed)

# should we modify the object and remove the original counts?

# write anndata to h5
scanorama_integrated_adata.write(filename = args.output_anndata)
