#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate snakemake-cia
mkdir -p workspace/
snakemake --directory workspace/ --configfile config/config.yaml -s CIAassembly_pipeline -j 16 all