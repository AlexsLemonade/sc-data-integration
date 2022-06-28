#!/bin/bash
set -euo pipefail

###############################################################################
# Preprocess the unfiltered SCE objects including, removing empty droplets and 
# creating the metadata file needed for downstream analyses. Then run 
# scpca-downstream-analyses

# **Note: Before running this, the original loom files must be converted to SCE 
# objects by running `00-convert-loom.R`.

# Usage, note that the --downstream_repo option is required as there is no default set
# bash run-downstream-analyses.sh \
#  --downstream_repo "full path to scpca-downstream-analyses repo"

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
#   scpca-downstream-analyses will be used.

###############################################################################

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

processed_library_df=$script_directory/../sample-info/hca-processed-libraries.tsv
metadata_file=$script_directory/../sample-info/hca-library-metadata.tsv
unfiltered_sce_dir=$script_directory/../data/human_cell_atlas/sce
filtered_sce_dir=$script_directory/../results/human_cell_atlas/filtered_sce
downstream_metadata_file=$script_directory/../sample-info/hca-downstream-metadata.tsv
results_dir=$script_directory/../results/human_cell_atlas/scpca-downstream-analyses
mito_file=$script_directory/../reference-files/gencode.v27.mitogenes.txt

# grab variables from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# run Rscript to generate metadata file 
Rscript --vanilla 01-preprocess-sce.R \
  --library_file $processed_library_df \
  --unfiltered_sce_dir $unfiltered_sce_dir \
  --filtered_sce_dir $filtered_sce_dir \
  --output_metadata $downstream_metadata_file

# check for Snakefile in downstream repo 
if [ ! -f $downstream_repo/Snakefile ]; then
  echo "The path provided for `--downstream_repo` is missing a Snakefile. 
        Double check you have provided the correct path." 
  exit 1
fi

# activate snakemake environment before running snakemake 
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate snakemake
snakemake --cores 2 \
  -s $downstream_repo/Snakefile \
  --config results_dir=$results_dir \
  project_metadata=$downstream_metadata_file \
  mito_file=$mito_file
  