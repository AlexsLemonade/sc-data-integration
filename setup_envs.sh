#!/bin/bash
#
# Sets up conda environments for snakemake
# On Apple Silicon devices, uses Intel R for Bioconductor compatibility

set -euo pipefail

# Run from the script file location
cd "$(dirname "${BASH_SOURCE[0]}")"

# If on OSX with Apple Silicon, build the Intel version of R
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  echo "On Apple Silicon: Building R environment and required R packages with CONDA_SUBDIR=osx-64 for compatibility"
  CONDA_SUBDIR=osx-64 snakemake --use-conda --conda-create-envs-only -c1 build_renv
fi

echo "Building conda environments"
snakemake --use-conda --conda-create-envs-only -c1
