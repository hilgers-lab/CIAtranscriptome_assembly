# CIAtranscriptome_assembly

Snakefile pipeline of all the steps taken to reconstruct the CIA transcriptome assembly as in (Alfonso-Gonzalez, 2022). This pipeline is based on Drosophila genome. 

## Requirements 

Packages dependencies:

* [FLAIR.V1.1](https://github.com/BrooksLabUCSC/flair/tree/v1.0)
* [SQANTI3](https://github.com/ConesaLab/SQANTI3)
* [minimap2](https://github.com/lh3/minimap2)

R dependencies: 
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html)  
* [tidyverse](https://tidyverse.tidyverse.org/) 
R Bioconductor: 
* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) 
* [pylranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html)
* [BSgenome](https://www.rdocumentation.org/packages/BSgenome/versions/1.40.1/topics/getSeq-methods) 


## Run 

Pipeline takes as input FASTQ files obtained from long read sequencing experiments in a folder named FASTQ

```
snakemake -s CIAassembly_pipeline  
```

