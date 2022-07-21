#!/usr/bin/env python3

import anndata as adata
import scanorama
import scipy

# input to function should be merged hdf5 object
# output should be anndata object with corrected gene expression and embeddings 

def integrate_scanorama(merged_adata):

    # grab variable gene list from merged object 
    var_genes = list(merged_adata.uns['variable_genes'])
    
    # subset merged object to only contain variable genes
    merged_adata = merged_adata[merged_adata.obs_names, var_genes]
    
    # create a dictionary with sample and associated cells for that sample
    sample_dict = merged_adata.obs.groupby("batch").indices
    sample_names = sample_dict.keys()
    
    # split merged object into a list of matrices, the logcounts subset to HVG  
    split_logcounts = []
    split_genes = []
    split_adata = []
    for sample in sample_names:
        anndata = merged_adata[sample_dict[sample]]
        split_adata.append(anndata)
        split_genes.append(anndata.var_names)
        
        logcounts = anndata.layers["logcounts"]
        logcounts = scipy.sparse.csr_matrix(logcounts)
        split_logcounts.append(logcounts)
    
    # perform integration with returning embeddings (SVD)
    integrated, corrected, genes = scanorama.correct(split_logcounts, 
                                                     split_genes,
                                                     return_dimred = True)
    
    # add corrected gene expression to individual anndata objects 
    for idx, anndata in enumerate(split_adata): 
        anndata.obsm["scanpy_SVD"] = integrated[idx]
        anndata.layers["scanpy_corrected"] = corrected[idx]
        
    # merge anndata back together 
    integrated_anndata_obj = adata.concat(split_adata, merge = "same")
    
    return integrated_anndata_obj
    



