#!/bin/bash
set -euo pipefail

###############################################################################
# Preprocess the unfiltered SCE objects including, removing empty droplets and
# creating the metadata file needed for downstream analyses. Then run
# scpca-downstream-analyses

# **Note: Before running this, the original loom files must be converted to SCE
# objects by running `00-obtain-sce.R`.

# Usage, note that the --downstream_repo option is required as there is no default set
# bash 01-run-downstream-analyses.sh \
#  --downstream_repo "path to scpca-downstream-analyses repo"

# Other parameters include:
# --processed_library_df: path to the file listing all libraries that should be
#   included in processing.
# --unfiltered_sce_dir: path to folder where all unfiltered sce objects are located
# --filtered_sce_dir: path to folder where all filtered sce objects are to be stored
# --downstream_metadata_file: Path to write metadata file to be used to run
#   scpca-downstream-analyses
# --results_dir: Path to save results from running scpca-downstream-analyses
# --mito_file: Path to file with list of mitochondrial genes to use.
#   If not specified, the mitochondrial file for Ensembl-104 present in
#   scpca-downstream-analyses will be used
# --repeat_filtering: An option to repeat the filtering of empty droplets if filtering files already exist.
#   Default here is not to repeat filtering. To turn on use `--repeat_filtering yes`
# --s3_bucket: S3 bucket to sync output from scpca-downstream-analyses.
#   Syncing will only be performed if a bucket is provided
# --cores: Default number of CPU cores to use for snakemake

###############################################################################


call_dir="$(pwd)"

# included paths are relative to the project root, which we find relative to this script
script_dir="$(dirname "${BASH_SOURCE[0]}")"
cd "$script_dir" || exit

# grab full path to project root
project_root=$(git rev-parse --show-toplevel)

# change back to calling directory for any command line paths
cd "$call_dir" || exit

processed_library_df=${project_root}/sample-info/hca-processed-libraries.tsv
metadata_file=${project_root}/sample-info/hca-library-metadata.tsv
unfiltered_sce_dir=${project_root}/data/human_cell_atlas/sce
filtered_sce_dir=${project_root}/results/human_cell_atlas/filtered_sce
downstream_metadata_file=${project_root}/sample-info/hca-downstream-metadata.tsv
results_dir=${project_root}/results/human_cell_atlas/scpca-downstream-analyses
mito_file=${project_root}/reference-files/gencode.v27.mitogenes.txt
repeat_filtering="no"
s3_bucket=""
cores=2

# grab variables from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# create flag for repeat_filtering
if [[ $repeat_filtering == "yes" ]]; then
  repeat_filtering_flag="--repeat_filtering"
else
  repeat_filtering_flag=""
fi

# run Rscript to generate metadata file
Rscript --vanilla ${script_dir}/utils/preprocess-sce.R \
  --library_file $processed_library_df \
  --unfiltered_sce_dir $unfiltered_sce_dir \
  --filtered_sce_dir $filtered_sce_dir \
  --output_metadata $downstream_metadata_file \
  $repeat_filtering_flag

# check for Snakefile in downstream repo
if [[ ! -f $downstream_repo/Snakefile ]]; then
  echo "The path provided for '--downstream_repo' is missing a Snakefile.
        Double check you have provided the correct path."
  exit 1
fi

# run downstream analysis workflow from scpca-downstream-analyses
cd $downstream_repo

# check if Apple Silicon, and build environments if required
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  CONDA_SUBDIR=osx-64 snakemake --cores $cores \
  --use-conda --conda-create-envs-only \
  build_renv
fi

snakemake --cores $cores \
  --use-conda \
  --config results_dir=$results_dir \
  project_metadata=$downstream_metadata_file \
  mito_file=$mito_file
cd $call_dir


# sync output from snakefile to aws
if [[ -n $s3_bucket ]]; then
  aws s3 sync $results_dir $s3_bucket
fi
