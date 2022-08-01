#!/usr/bin/env python3

import anndata as adata
import scanorama
import scipy
import sys
import random


def integrate_scanorama(merged_adata,
                        batch_column='batch',
                        seed=None):
    """
    Function to integrate single cell gene expression data stored in an AnnData
    object. The AnnData object should contain the data from all libraries to be
    integrated and contain a column in the `anndata.obs` that indicates
    which library or experiment each cell is originally from. This can be specified
    using the `batch_column` argument.


    Parameters
    ----------
    merged_adata : AnnData object that contains the merged data from all libraries
        to be integrated. This object must include a list of highly variable genes
        found in `anndata.uns['variable_genes']`.
    batch_column : The name of the column in `anndata.obs` that indicates the batches
        for each cell, typically this corresponds to the library id.
    seed : Random seed to set prior to scanorama. A seed will only be set if this
        is not `None`.

    Returns
    -------
    integrated_anndata_obj : AnnData object containing the highly variable genes and
        `scanorama_SVD` in the `anndata.obsm` and `scanorama_corrected` gene
        expression matrix

    """

    # set seed
    random.seed(seed)

    # grab variable gene list from merged object
    try:
        var_genes = list(merged_adata.uns['variable_genes'])
    except KeyError:
        print("Variable genes cannot be found in anndata object."
              " Make sure they are stored in adata.uns['variable_genes'].",
              file = sys.stderr)
        raise

    # subset merged object to only contain variable genes
    merged_adata = merged_adata[merged_adata.obs_names, var_genes]

    # create a dictionary with sample and associated cells for that sample
    try:
        library_dict = merged_adata.obs.groupby(batch_column).indices
    except KeyError:
        print("Provided batch_column cannot be found in anndata object."
              f" Make sure it is stored in adata.obs[{batch_column}].",
              file = sys.stderr)
        raise

    # split merged object into a list of matrices corresponding to
    # the logcounts subset to HVG
    split_logcounts = []
    split_genes = [] # need a list of gene lists for input to scanorama
    split_adata = [] # also need to split the anndata object to add the integrated results back
    for library in library_dict:
        # subset merged anndata with cells for each library
        anndata = merged_adata[library_dict[library]]
        split_adata.append(anndata)

        # add to list of genes
        split_genes.append(anndata.var_names)

        # split out logcounts and convert to required csr_matrix
        logcounts = anndata.layers["logcounts"]
        logcounts = scipy.sparse.csr_matrix(logcounts)
        split_logcounts.append(logcounts)

    # perform integration with returning embeddings (SVD)
    integrated, corrected, genes = scanorama.correct(split_logcounts,
                                                     split_genes,
                                                     return_dimred = True)

    # add corrected gene expression and embeddings to individual anndata objects
    for idx, anndata in enumerate(split_adata):
        anndata.obsm["scanorama_SVD"] = integrated[idx]
        anndata.layers["scanorama_corrected"] = corrected[idx]

    # merge anndata back together into one integrated object with original and integrated data
    integrated_anndata_obj = adata.concat(split_adata, merge = "same")

    # add unstructured metadata from original merged object back to integrated object
    integrated_anndata_obj.uns = merged_adata.uns

    return integrated_anndata_obj
