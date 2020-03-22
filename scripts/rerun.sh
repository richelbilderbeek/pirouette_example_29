#!/bin/bash
#
# Re-run the code locally, to re-create the data and figure.
#
# Usage:
#
#   ./scripts/rerun.sh
#
#SBATCH --partition=gelifes
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --job-name=pirex29
#SBATCH --output=example_29.log
#
rm -rf example_29
rm *.png
time Rscript example_29.R
zip -r pirouette_example_29.zip example_29 example_29.R scripts *.png

