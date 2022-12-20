#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate snakemake-cia
mkdir -p workspace/
snakemake --use-conda --directory workspace/ --configfile config/config.yaml -s CIAassembly_pipeline -j 16 all
#snakemake --profile ~/.config/mpi-ie-slurm --directory workspace/ --configfile config/config-mpi.yaml -s CIAassembly_pipeline -j 16 all
