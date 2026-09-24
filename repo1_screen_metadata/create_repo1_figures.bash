#!/bin/bash

# Build the REPO1 compound metadata tables and create the REPO1 compound figures.
#
#   1. raw_data/   the metadata notebooks (Python) turn REPO1_April2025.xlsx into tidy annotation tables:
#                  the mechanism of action (MOA) annotation, the target annotation with Entrez IDs, and the genomic
#                  positions of the target genes. The Entrez ID and position lookups use the MyGene.info web service
#                  (only gene symbols / Entrez IDs are sent), so this step needs an internet connection.
#   2. notebooks/  the figure notebooks (R) create the figures (png and pdf) in notebooks/figures/.
#
# Requirements
#   step 1: python3 with jupyter (nbconvert), pandas, openpyxl and mygene
#   step 2: R with dplyr, tidyr, stringr, readxl, forcats, scales, ggplot2 and gganatogram (not on CRAN or conda-forge;
#           install it with devtools::install_github("jespermaag/gganatogram")), and jupyter (nbconvert)
#
# Usage:
#   bash create_repo1_figures.bash
# Optional environment variables:
#   PY_ENV, R_ENV     conda environments to activate for step 1 and step 2 (default: use the current environment)
#   SKIP_METADATA=1   skip step 1 and use the tables already in raw_data/ (no internet connection needed), e.g.
#                     SKIP_METADATA=1 bash create_repo1_figures.bash

# stop at the first error
set -e

here="$(cd "$(dirname "$0")" && pwd)"

# activate a conda environment, if one is named
activate_env() {
    if [ -n "$1" ]; then
        eval "$(conda shell.bash hook)"
        conda activate "$1"
    fi
}

# ---- step 1: metadata processing ----
if [ "${SKIP_METADATA:-0}" != "1" ]; then
    activate_env "${PY_ENV:-}"
    cd "${here}/raw_data"

    # the order matters: the gene positions are looked up for the genes of the target annotation
    for notebook in build_repo1_moa_annotation build_repo1_target_annotation build_repo1_target_gene_positions; do
        # convert the notebook to a script in the nbconverted folder and run it (the notebooks use paths relative to raw_data/)
        jupyter nbconvert --to script --output-dir=nbconverted/ "${notebook}.ipynb"
        echo "Running nbconverted/${notebook}.py"
        python3 "nbconverted/${notebook}.py"
    done
fi

# ---- step 2: figures ----
activate_env "${R_ENV:-}"
cd "${here}/notebooks"   # the notebooks use paths relative to their own folder

# convert all notebooks to script files into the nbconverted folder
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run each script to create its figures in the figures folder
for script in nbconverted/*.r; do
    echo "Running ${script}"
    Rscript "${script}"
done

# drawing the plots to the screen device writes Rplots.pdf, which is not needed
rm -f Rplots.pdf

echo "Figures are in $(pwd)/figures"
