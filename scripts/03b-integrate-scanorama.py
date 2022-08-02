#!/usr/bin/env python3
"""
Script to perform integration on a merged AnnData object using Scanorama

Takes an HDF5 file as input containing a merged AnnData object across at least 2 single-cell libraries
and performs correction with scanorama. The output is an HDF5 file containing the integrated AnnData object.
Requires a `batch_column` to be present as a column in the `obs` of the AnnData which corresponds to the
original library ID for each cell. To return only the corrected data, removing the `X` matrix containing
the raw counts and the `logcounts` layer containing the log-normalized data, use the `--corrected_only` flag.

"""


import os
import anndata as adata
import argparse
import re
from utils.integrate_scanorama import integrate_scanorama

# define arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--input_anndata',
                    dest = 'input_anndata',
                    help = 'Path to HDF5 file with merged AnnData object to integrate')
parser.add_argument('-i', '--output_anndata',
                    dest = 'output_anndata',
                    help = 'Path to HDF5 file to save the integrated AnnData object')
parser.add_argument('-b', '--batch_column',
                    dest = 'batch_column',
                    default = 'batch',
                    help = ('The name of the column in `anndata.obs` that indicates the batches for each cell, '
                            ' typically this corresponds to the library id.'))
parser.add_argument('--corrected_only',
                    dest = 'corrected_only',
                    action = 'store_true',
                    help = 'Boolean indicating to return only the corrected data or all raw data.'
                    ' Default will return all data. To return only corrected data, use --corrected_only.')
parser.add_argument('-s', '--seed',
                    dest = 'seed',
                    default = None,
                    help = 'Random seed to set for scanorama.')

args = parser.parse_args()

# compile extension regex
file_ext = re.compile(r"\.hdf5$|.h5$", re.IGNORECASE)

# check that input file exists, if it does exist, make sure it's an h5 file
if not os.path.exists(args.input_anndata):
    raise FileExistsError("--input_anndata file not found.")
elif not file_ext.search(args.input_anndata):
    raise ValueError("--input_anndata must end in either .hdf5 or .h5 and contain a merged AnnData object.")

# make sure output file path is h5 file
if not file_ext.search(args.output_anndata):
    raise ValueError("--output_anndata must provide a file path ending in either .hdf5 or .h5.")

# check that output file directory exists and create directory if doesn't exist
integrated_adata_dir = os.path.dirname(args.output_anndata)
if not os.path.isdir(integrated_adata_dir):
    os.makedirs(integrated_adata_dir, exist_ok = True)

# read in merged anndata object
merged_adata = adata.read_h5ad(args.input_anndata)

# integrate anndata with scanorama
scanorama_integrated_adata = integrate_scanorama(merged_adata,
                                                 batch_column = args.batch_column,
                                                 seed = args.seed)

# remove raw data and logcounts to minimize space if corrected_only is true
if args.corrected_only:
    print("Removing raw data and log-normalized data. Only corrected data will be returned.")
    scanorama_integrated_adata.X = None
    del scanorama_integrated_adata.layers["logcounts"]

# write anndata to h5
scanorama_integrated_adata.write(filename = args.output_anndata)
