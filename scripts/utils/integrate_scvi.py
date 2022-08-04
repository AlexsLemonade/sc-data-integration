#!/usr/bin/env python3

import sys
import scvi

def integrate_scvi(merged_adata,
                   integrate_genes,
                   batch_column = ['batch'],
                   categorical_covariate_columns = None,
                   continuous_covariate_columns = ['subsets_mito_percent'],
                   seed=None):

    """
    Function to integrate single cell gene expression data stored in an AnnData
    object using `scVI`. The AnnData object should contain the data from all libraries to be
    integrated and contain a column in the `anndata.obs` that indicates
    which library or experiment each cell is originally from. This can be specified
    using the `batch_column` argument. Can also indicate additional covariates to consider,
    such as patient, gender, etc using the `categorical_covariate_columns` argument. Continuous covariates,
    such as mitochondrial percentage are used by default. If different continuous covariates
    are to be used, like ribosomal percentage, they can be indicated using the
    `continuous_covariate_columns` argument.


    Parameters
    ----------
    merged_adata : AnnData object that contains the merged data from all libraries
        to be integrated. This object must include a list of highly variable genes
        found in `anndata.uns['variable_genes']`.
    integrate_genes : List of genes to use for integration. AnnData objects will be
        subset to contain only these genes before integration. Typically this corresponds
        to a list of highly variable genes.
    batch_column : The name of the column in `anndata.obs` that indicates the batches
        for each cell, typically this corresponds to the library id. Default is `batch`.
    categorical_covariate_columns : A list of columns containing additional categorical
        data to be included as a covariate. Default is None.
    continuous_covariate_columns : A list of columns containing additional continous
        data to be included as a covariate. Default is `['subsets_mito_percent']`.
    seed : Random seed to set prior to scanorama. A seed will only be set if this
        is not `None`.

    Returns
    -------
    integrated_anndata_obj : AnnData object containing the highly variable genes and
        `scvi_latent` in the `anndata.obsm` and `scvi_corrected` gene
        expression matrix

    """
    # set seed
    scvi.settings.seed = seed

    # make sure type of integrate_genes is a list
    if not type(integrate_genes) is list:
        raise TypeError("`integrate_genes` must be a list.")

    # subset merged object to only contain genes used for integration
    try:
        merged_adata = merged_adata[merged_adata.obs_names, integrate_genes]
    except KeyError:
        print("Genes provided in `integrate_genes` are not present as rows in the AnnData object.",
              file = sys.stderr)
        raise


    # if additional covariate columns are provided concatentate with batch column
    if type(batch_column) != list:
        batch_column = [batch_column]
    if len(batch_column) != 1:
        raise ValueError("Please provide only one column for batch_column."
                         " If any other covariates should be included, include them as categorical_covariate_columns.")


    if categorical_covariate_columns is None:
        categorical_covariate_columns = batch_column
    elif type(categorical_covariate_columns) == str:
        # if only one covariate is provided, need to convert to list first
        categorical_covariate_columns = batch_column + [categorical_covariate_columns]
    elif type(categorical_covariate_columns) == list:
        categorical_covariate_columns = batch_column + categorical_covariate_columns
    else:
        raise TypeError("categorial_covariate_columns must be either type str or list.")

    # make sure that input for covariate columns is a list
    if type(continuous_covariate_columns) != list:
        continuous_covariate_columns = [continuous_covariate_columns]

    # check that batch, categorical, and continuous categorical keys are all found in obs
    # batch is now stored in categorical covariates so don't need to check separately
    missing_cols = [col for col in (continuous_covariate_columns + categorical_covariate_columns)
                    if col not in merged_adata.obs]
    if missing_cols:
        raise ValueError("The following columns from batch_column, categorical_covariate_columns,"
                         " and/or continuous_covariate_columns are not present in adata.obs:",
                         ", ".join(missing_cols))

    # make sure that adata is present as a copy and not just a view
    merged_adata = merged_adata.copy()
    # convert raw counts to csr matrix
    merged_adata.X = merged_adata.X.tocsr()

    # create an scvi model object
    # raw counts are required as input
    scvi.model.SCVI.setup_anndata(merged_adata,
                                  categorical_covariate_keys = categorical_covariate_columns,
                                  continuous_covariate_keys = continuous_covariate_columns)

    # create and train model
    scvi_model = scvi.model.SCVI(merged_adata)
    scvi_model.train()

    # dimensionality reduction, mean of the approximate posterior
    # save as latent, but recommended to treat like PCA and use as input to UMAP
    merged_adata.obsm["scvi_latent"] = scvi_model.get_latent_representation()

    # normalized/corrected gene expression data
    merged_adata.layers["scvi_corrected"] = scvi_model.get_normalized_expression()

    return merged_adata

