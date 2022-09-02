# Analysis Templates

## Single group integration comparison template

The `01-single-group-integration-check-template.Rmd` contains a template `Rmd` file that can be used to evaluate a single integration method on a single group of libraries.
The template requires a pre-integrated and a post-integrated `SingleCellExperiment` object with all libraries that were integrated contained in the same object.
The object pre-integrated object must contain principal components calculated prior to integration labeled as `PCA` along with UMAP embeddings labeled as `UMAP`, and the post-integrated object must contain principal components or PC-equivalent embeddings labeled as `<integration_method>_PCA` along with UMAP embeddings labeled as `<integration_method>_UMAP`.
Both objects must contain a column in the `colData` indicating which libary each cell originated from, typically labeled `batch`, and a column containing the celltype information.

The following are a list of parameters in the template notebook that can be used to evaluate integration:

- `integration_method`: Which method was used for integration, one of `fastmnn`, `harmony`, `scvi`, `scanorama`, `cca`, or `rpca`
- `batch_column`: The column in the `colData` used to indicate which library each cell originated from, default is `batch`
- `celltype_column`: The column in the `colData` used to indicate which cell type each cell belongs to, default is `celltype`
- `num_pcs`: Number of principal components to use as input for calculating kBET, batch silhouette width, and batch ARI, default is 20
- `ari_k_min`: For calculating the batch ARI, min number of centers to use for clustering by k-means, default is 5
- `ari_k_max`: For calculating the batch ARI, max number of centers to use for clustering by k-means, default is 25
- `ari_k_increment`: For calculating the batch ARI, increment to increase the number of centers to use for clustering by k-means, e.g. for the default `ari_k_increment` of 5 with default `ari_k_min` and `ari_k_max`, a sequence will be created of `c(5, 10, 15, 20, 25)`
- `k0_fraction_min`: For calculating kBET, min fraction of sample size to set k0, default is 0.05
- `k0_fraction_max`: For calculating kBET, max fraction of sample size to set k0, default is 0.25
- `k0_fraction_increment`: For calculating kBET, increment to increase the fraction of sample size to use for setting k0, e.g. for the default `k0_fraction_increment` of 0.05 with default `k0_fraction_min` and `k0_fraction_max`, a sequence will be created of `c(0.05, 0.01, 0.15, 0.20, 0.25)`
- `num_regression_pcs`: For performing PCA regression for batch effect, the number of PCs to use, default is 50
- `sig_threshold`: For performing PCA regression for batch effect, threshold for considering a corrected P-value from PC~batch regression as significant, default is 0.05
