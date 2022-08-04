#!/usr/bin/env python3
"""
Script to perform integration on a merged AnnData object using scVI

Takes an HDF5 file as input containing a merged AnnData object across at least 2 single-cell libraries
and performs correction with scVI. The output is an HDF5 file containing the integrated AnnData object.
Requires a `batch_column` to be present as a column in the `obs` of the AnnData which corresponds to the
original library ID for each cell. Can also provide either `continuous_covariates` or `categorical_covariates`
to be considered for integration. The provided covariates must be included as columns in the `obs` of the AnnData.
To return only the corrected data, removing the `X` matrix containing
the raw counts and the `logcounts` layer containing the log-normalized data, use the `--corrected_only` flag.
To use only highly variable genes for integration, use the `--use_hvg` flag. When this flag is used the returned object
will only contain highly variable genes.
"""

import os
import anndata as adata
import argparse
import re
from utils.integrate_scvi import integrate_scvi

# define arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_anndata',
                    dest = 'input_anndata',
                    required = True,
                    help = 'Path to HDF5 file with merged AnnData object to integrate')
parser.add_argument('-o', '--output_anndata',
                    dest = 'output_anndata',
                    required = True,
                    help = 'Path to HDF5 file to save the integrated AnnData object')
parser.add_argument('-b', '--batch_column',
                    dest = 'batch_column',
                    default = 'batch',
                    help = ('The name of the column in `anndata.obs` that indicates the batches for each cell, '
                            ' typically this corresponds to the library id.'))
parser.add_argument('--categorical_covariates',
                    dest = 'categorical_covariates',
                    default = None,
                    help = 'A comma-separated list of columns containing additional categorical data to be'
                    ' included as a covariate. Default is None.')
parser.add_argument('--continuous_covariates',
                    dest = 'continuous_covariates',
                    default = "subsets_mito_percent",
                    help = 'A comma-separated list of columns containing additional continous data to be'
                    ' included as a covariate. Default is "subsets_mito_percent".')
parser.add_argument('--use_hvg',
                    dest = 'use_hvg',
                    action = 'store_true',
                    default = False,
                    help = 'Boolean indicating whether or not to use only highly variable genes for data integration.'
                    'If --use_hvg is used, the returned integrated object will only contain the highly variable genes.')
parser.add_argument('--corrected_only',
                    dest = 'corrected_only',
                    action = 'store_true',
                    help = 'Boolean indicating to return only the corrected data or all raw data.'
                    ' Default will return all data. To return only corrected data, use --corrected_only.')
parser.add_argument('-s', '--seed',
                    dest = 'seed',
                    type=int,
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
os.makedirs(integrated_adata_dir, exist_ok = True)

# read in merged anndata object
merged_adata = adata.read_h5ad(args.input_anndata)

if args.use_hvg:
    try:
        var_genes = list(merged_adata.uns['variable_genes'])
    except KeyError:
        print("Variable genes cannot be found in anndata object."
              " If using --use_hvg, make sure HVG are stored in adata.uns['variable_genes'].")
        raise
else:
    var_genes = list(merged_adata.var_names)

# split covariates from comma separated strings into lists and check type
if not args.categorical_covariates is None:
    if type(args.categorical_covariates) is str:
        categorical_covariates = [covariate.strip() for covariate in args.categorical_covariates.split(',')]
    else:
        raise TypeError("--categorical_covariates must be a comma separated list")
else:
    categorical_covariates = None


if not args.continuous_covariates is None:
    if type(args.continuous_covariates) is str:
        continuous_covariates = [covariate.strip() for covariate in args.continuous_covariates.split(',')]
    else:
        raise TypeError("--continuous_covariates must be a comma separated list")
else:
    continuous_covariates = None

# integrate anndata with scvi
scvi_integrated_adata = integrate_scvi(merged_adata,
                                       integrate_genes = var_genes,
                                       batch_column = args.batch_column,
                                       categorical_covariate_columns = categorical_covariates,
                                       continuous_covariate_columns = continuous_covariates,
                                       seed = args.seed)

# remove raw data and logcounts to minimize space if corrected_only is true
if args.corrected_only:
    print("Removing raw data and log-normalized data. Only corrected data will be returned.")
    scvi_integrated_adata.X = None
    del scvi_integrated_adata.layers["logcounts"]


# write anndata to h5
scvi_integrated_adata.write(filename = args.output_anndata)
