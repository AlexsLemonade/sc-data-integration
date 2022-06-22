#!/bin/bash
set -euo pipefail

###############################################################################
# Create the metadata file needed for downstream analyses and then run 
# scpca-downstream-analyses

# Usage, note that the --snakefile option is required as there is no default set
# bash run-downstream-analyses.sh \
#  --downstream_repo "full path to scpca-downstream-analyses repo"

# Other parameters include: 
# --processed_library_df: The path to the file listing all libraries that should be included 
#   in the conversion from loom to SCE. This file must contain the 
#   `library_biomaterial_id` column. 
# --metadata_file: The path to the metadata file for all libraries. This file must 
#   contain columns for `library_biomaterial_id` and `sample_biomaterial_id`
# --sce_dir: Path to the folder where all SCE objects are saved locally
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
sce_dir=$script_directory/../data/human_cell_atlas/sce
downstream_metadata_file=$script_directory/../sample-info/hca-downstream-metadata.tsv
results_dir=$script_directory/../results/scpca-downstream-analyses

# grab variables from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# use mito file provided in scpca-downstream-analyses if not provided
mito_file=${mito_file:-"${downstream_repo}/reference-files/Homo_sapiens.GRCh38.104.mitogenes.txt"}

# run Rscript to generate metadata file 
Rscript --vanilla 01-preprocess-sce.R \
  --library_file $processed_library_df \
  --full_metadata_file $metadata_file \
  --sce_dir $sce_dir \
  --output_metadata $downstream_metadata_file

# paths in the downstream metadata file are relative to the sce directory 
# move to sce directory before running the snakefile 
cd $sce_dir
mkdir -p $results_dir

# activate snakemake environment before running snakemake 
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate snakemake
snakemake --cores 2 \
  -s $downstream_repo/Snakefile \
  --config results_dir=$results_dir \
  project_metadata=$downstream_metadata_file \
  mito_file=$mito_file
  