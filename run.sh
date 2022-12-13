#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate snakemake
mkdir -p workspace/
snakemake --profile ~/.config/mpi-ie-slurm --directory workspace/ --configfile config/config-mpi.yaml -s CIAassembly_pipeline -j 16 all