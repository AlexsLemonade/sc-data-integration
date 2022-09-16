# Analysis Templates

## Single group integration comparison template

The `01-single-group-integration-check-template.Rmd` contains a template `Rmd` file that can be used to evaluate all tested integration methods on a single group of libraries.
The template reads in the directories that hold the pre-integrated and post-integrated `SingleCellExperiment` objects with all libraries that were integrated contained in each object.
The template also requires a `group_name` corresponding to the name used to label the RDS files containing the integrated `SingleCellExperiment` objects.
The pre-integrated object must be located within the directory provided for the pre-integrated object and have the following file pattern: `<group_name>_merged_sce.rds`.
The post-integrated objects must all be located within the directory provided for the post-integrated object and have both the `group_name` and one of the following required integration methods in the file name: `fastmnn`, `harmony`, `cca`, `rpca`, `scanorama`, and `scvi`.

The pre-integrated object must contain principal components calculated prior to integration labeled as `PCA` along with UMAP embeddings labeled as `UMAP`, and the post-integrated objects must contain principal components or PC-equivalent embeddings labeled as `<integration_method>_PCA` along with UMAP embeddings labeled as `<integration_method>_UMAP`.
Both objects must contain a column in the `colData` indicating which libary each cell originated from, typically labeled `batch`, and a column containing the celltype information.

The following are a list of parameters in the template notebook that can be used to evaluate integration:

- `merged_sce_dir`: The directory containing the pre-integrated `SingleCellExperiment` object stored as an RDS file and named with the `group_name`
- `integrated_sce_dir`: The directory containing all post-integrated `SingleCellExperiment` objects stored as RDS files named with both the `group_name` and one of the required names for `integration_method`
- `group_name`: The name used to label the group of libraries that make up the dataset
- `batch_column`: The column in the `colData` used to indicate which library each cell originated from, default is `batch`
- `celltype_column`: The column in the `colData` used to indicate which cell type each cell belongs to, default is `celltype`
- `num_pcs`: Number of principal components to use as input for calculating kBET, batch silhouette width, and batch ARI, default is 20
- `ari_k_min`: For calculating the batch ARI, min number of centers to use for clustering by k-means, default is 5
- `ari_k_max`: For calculating the batch ARI, max number of centers to use for clustering by k-means, default is 25
- `ari_k_increment`: For calculating the batch ARI, increment to increase the number of centers to use for clustering by k-means, e.g. for the default `ari_k_increment` of 5 with default `ari_k_min` and `ari_k_max`, a sequence will be created of `c(5, 10, 15, 20, 25)`
- `num_regression_pcs`: For performing PCA regression for batch effect, the number of PCs to use, default is 50
- `sig_threshold`: For performing PCA regression for batch effect, threshold for considering a corrected P-value from PC~batch regression as significant, default is 0.05

The following metrics will be calculated in the notebook:

- iLISI as defined in [Korsunksky et al.](https://doi.org/10.1038/s41592-019-0619-0)
- batch adjusted rand index (ARI)
- batch average silhouette width (ASW)
- PC batch variance and scaled PC regression, derived from the [manuscript introducing kBET](https://doi.org/10.1038/s41592-018-0254-1)

Note that we no longer include [kBET](https://github.com/theislab/kBET) as a metric because it is highly sensitive to bias present in the data.
Testing of kBET revealed that in almost all scenarios the observed rejection rate was close to 1 indicating batch effect was present in the data.
