#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate snakemake-cia
WORKSPACE=workspace-test
mkdir -p $WORKSPACE
#snakemake --directory workspace/ --configfile config/config.yaml -s CIAassembly_pipeline -j 16 all
snakemake --directory $WORKSPACE --configfile config/config-test.yaml -s CIAassembly_pipeline -j 8 all