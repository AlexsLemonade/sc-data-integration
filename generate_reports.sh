#!/bin/bash
#
# Generate reports for HCA and simulation data
# This script assumes conda environments have been built with `setup_envs.sh`, and that `snakemake` is in the PATH

set -euo pipefail


# Generate HCA project reports using the default config.yaml
snakemake -c1 --use-conda


# Generate simulation project reports using provided parameters
snakemake -c1 --use-conda \
  --config processed_tsv="sample-info/scib-simulated-processed-libraries.tsv" \
  sce_dir="data/scib_simulated/sce" \
  results_dir="results/scib_simulated" \
  add_celltype=false \
  pepfile="sample-info/scib-simulated-project-pep.yaml" \
  continuous_covariates=None

