#!/bin/sh

# Check if we are on an Apple Silicon mac
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  if [[ $CONDA_SUBDIR == 'osx-64' ]]; then
    # Record CONDA_SUBDIR if it is set
    conda env config vars set CONDA_SUBDIR=osx-64
  else
    # Error if CONDA_SUBDIR was not set
    echo 'BioConductor is not compatible with arm64, rerun with `CONDA_SUBDIR=osx-64`' >&2
    echo ' e.g. `CONDA_SUBDIR=osx-64 snakemake --use-conda --conda-create-envs-only -c1`' >&2
    exit 1
  fi
fi



# install packages required for renv
Rscript --no-init-file -e 'install.packages(c("jsonlite","purrr"), repos = "https://cloud.r-project.org")'
# restore packages
Rscript -e 'renv::restore()'
