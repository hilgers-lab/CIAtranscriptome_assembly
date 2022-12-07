#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate snakemake
snakemake --profile ~/.profile/mpi-ie-slurm --workdir workspace/ --configfile config/config.yaml -j 16 -s CIAassembly_pipeline --   