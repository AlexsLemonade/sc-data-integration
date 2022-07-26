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

    # set seed
    random.seed(seed)

    # need to input only HVG to build model so subset to HVG 
    # make sure variable gene list alread exists in merged object first
    try:
        var_genes = list(merged_adata.uns['variable_genes'])
        merged_adata = merged_adata[merged_adata.obs_names, var_genes]
    except KeyError:
        print("Variable genes cannot be found in anndata object."
              "Make sure they are stored in adata.uns['variable_genes'].")
        sys.exit(1)


    # if additional covariate columns are provided concatentate with batch column
    if type(batch_column) != list:
        batch_col = [batch_column]
    if len(batch_col) != 1:
        print("Please provide only one column for batch_column."
        "If any other covariates should be included, include them as categorical_covariate_columns.")
        sys.exit(1)


    if categorical_covariate_columns is not None:
        # if only one covariate is provided, need to convert to list first 
        if type(categorical_covariate_columns) == str:
            categorical_covariate_columns = batch_col + [categorical_covariate_columns]
        if type(categorical_covariate_columns) == list:
            categorical_covariate_columns = batch_col + categorical_covariate_columns
    else:
        categorical_covariate_columns = batch_col

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
    # if have additional covariates can add to `categorical_covariate_keys`
    scvi.model.SCVI.setup_anndata(merged_adata,
                                  categorical_covariate_keys = categorical_covariate_columns,
                                  continuous_covariate_keys = continuous_covariate_columns)

    # create and train model 
    scvi_model = scvi.model.SCVI(merged_adata)
    scvi_model.train()

    # dimensionality reduction, mean of the approximate posterior 
    merged_adata.obsm["scvi_PCA"] = scvi_model.get_latent_representation()

    # normalized/corrected gene expression data 
    merged_adata.layers["scvi_corrected"] = scvi_model.get_normalized_expression()

    return merged_adata
