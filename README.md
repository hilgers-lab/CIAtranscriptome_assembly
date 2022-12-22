# CIA Transcriptome Assembly

Snakefile pipeline of all the steps taken to reconstruct the CIA transcriptome assembly as in (Alfonso-Gonzalez, 2022). This pipeline is based on Drosophila genome. Input should be Drosophila based data. 

## Requirements 

All package dependencies are downloaded using conda, with the exception of SQANTI and FLAIR which are installed during the pipeline run using specific github commits. 

Before running, make sure you have conda installed, then run `conda create -n snakemake-cia -c conda-forge -c bioconda -c defaults snakemake` 
to create the pipeline environment. Then run snakemake with `--use-conda` to automatically create the environments during the pipeline run. 

Briefly, dependencies are listed here. 

### Packages dependencies:

* [FLAIR.V1.1](https://github.com/BrooksLabUCSC/flair/tree/v1.0)
* [SQANTI3](https://github.com/ConesaLab/SQANTI3)
* [minimap2](https://github.com/lh3/minimap2)

### R dependencies: 

* [optparse](https://cran.r-project.org/web/packages/optparse/index.html)  
* [tidyverse](https://tidyverse.tidyverse.org/) 

### R Bioconductor: 

* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) 
* [plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html)
* [BSgenome](https://www.rdocumentation.org/packages/BSgenome/versions/1.40.1/topics/getSeq-methods) 


## Before you run

Edit `config/config.yaml` to reflect the parameters you would like to use to run the pipeline, as well as
`config/units.tsv` to specify the sample `path` and `sample_type` -- one of `flam-seq`, `iso-seq`, `ont-cdna`, or `ont-direct`.
Data files can be gzipped or raw FASTA or FASTQ files. 

## Run

Modify the snakemake command in `run.sh` to use parameters that are appropriate
for your computing or cluster environment. The pipeline uses conda, so be sure to 
include `--use-conda` in the snakemake command. 

Execute `./run.sh`. 

## Testing

You can download a test dataset, for which the `config/units.tsv` is already configured, from Zenodo [here](https://doi.org/10.5281/zenodo.7438383). Run `tar -xzvf test.tar.gz` in this directory and run the pipeline using `./run.sh`. 

## Troubleshooting

Occasionally there can be an issue installing R packages in the `cia-sqanti`
environment. This will manifest in an error like this: 

```
Error in if (nzchar(SHLIB_LIBADD)) SHLIB_LIBADD else character() : 
argument is of length zero
```

if you run into this error, follow the instructions from [this thread](https://stackoverflow.com/questions/53813323/installing-r-packages-in-macos-mojave-error-in-if-nzcharshlib-libadd/54778735?stw=2#54778735).