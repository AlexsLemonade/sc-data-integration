#!/usr/bin/env python3

import sys
import random
import scvi
import scipy

seed=2022
random.seed(seed)

def integrate_scvi(merged_adata,
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

    # need to input only HVG to build model so subset to HVG 
    # make sure variable gene list alread exists in merged object first
    try:
        var_genes = list(merged_adata.uns['variable_genes'])
        merged_adata = merged_adata[merged_adata.obs_names, var_genes]
    except KeyError:
        print("Variable genes cannot be found in anndata object."
              " Make sure they are stored in adata.uns['variable_genes'].")
        raise


    # if additional covariate columns are provided concatentate with batch column
    if type(batch_column) != list:
        batch_column = [batch_column]
    if len(batch_column) != 1:
        raise ValueError("Please provide only one column for batch_column."
                         " If any other covariates should be included, include them as categorical_covariate_columns.")


    if categorical_covariate_columns is not None:
        # if only one covariate is provided, need to convert to list first 
        if type(categorical_covariate_columns) == str:
            categorical_covariate_columns = batch_column + [categorical_covariate_columns]
        if type(categorical_covariate_columns) == list:
            categorical_covariate_columns = batch_column + categorical_covariate_columns
    else:
        categorical_covariate_columns = batch_column

    # make sure that input for covariate columns is a list
    if type(continuous_covariate_columns) != list:
        continuous_covariate_columns = [continuous_covariate_columns]

    # check that batch, categorical, and continuous categorical keys are all found in obs
    # batch is now stored in categorical covariates so don't need to check separately
    obs_keys = set(merged_adata.obs.keys().tolist())
    if not obs_keys.intersection(set(categorical_covariate_columns)) or \
        not obs_keys.intersection(set(continuous_covariate_columns)):
        print("Make sure that batch_column, categorical_covariate_columns, and continuous_covariate_columns"
              " are all present as columns in adata.obs.")
        sys.exit(1)

    # convert raw counts to csr matrix 
    merged_adata.X = scipy.sparse.csr_matrix(merged_adata.X)
    # make sure that adata is present as a copy and not just a view
    merged_adata = merged_adata.copy()
    
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
