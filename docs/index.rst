.. CIA Transcriptome Assembly documentation master file, created by
   sphinx-quickstart on Wed Feb 22 16:03:05 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CIA Transcriptome Assembly's documentation!
======================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. image:: _static/dag.png

Overview
--------

The CIA Transcriptome Assembly pipeline is a custom Snakemake pipeline designed
to map long-read mRNA or cDNA-seq data and build transcriptome annotations from
in a robust manner. It accepts FASTQ files from ONT Direct RNA-seq, 
ONT PCR-cDNA, PacBio Iso-Seq, or FLAM-seq as input. 

It uses `FLAIR <https://github.com/BrooksLabUCSC/flair>`_, 
`SQANTI <https://github.com/ConesaLab/SQANTI>`_, and a custom end correction
step to build the annotations with high resolution.  

This software accompanies Alfonso-Gonzalez et al., Cell, 2023 *TODO INSERT
DOI HERE*. 

Software/hardware requirements
------------------------------

The CIA Transcriptome Assembly has been tested on CentOS 7.7 and Ubuntu 22.04 
LTS and should work on most Linux64 platforms. At least 32GB of memory are 
required for the STAR alignment step. 

Getting started
---------------

Download and install ``conda`` from 
`here <https://docs.conda.io/en/latest/miniconda.html>`_. Execute 

:: 

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    conda install mamba -n base -c conda-forge

Then execute
``make env`` from the base directory to make the base snakemake environment. 
This will result in an environment called ``snakemake-cia``. If you need
to remove the conda environment, execute ``make clean-env``. 

All other software will be installed via conda or during the runtime of the 
pipeline. 

Testing the pipeline
--------------------

To ensure that the pipeline can successfully run on your machine, we have 
provided a dataset consisting of four files, one from each experimental type,
that can be downloaded from `Zenodo <https://doi.org/10.5281/zenodo.7438383>`_.

To automatically download this dataset and run it, simply run ``make test`` 
which will run ``make env`` as described above if necessary, 
download and build the Drosophila annotation, the Zenodo dataset, 
and unpack and run it using the ``test.sh`` snakemake command. 

Configuration
-------------

Two files in the config directory control most of the behavior of the 
pipeline: ``units.tsv`` and ``config/config.yaml``. ``units.tsv`` contains
three columns: ``sample``, ``path``, and ``sample_type``. ``sample`` is simply
a unique identifier for the sample during the pipeline, ``path`` is 
the path to the demultiplexed FASTQ file of the respective sample (can be 
outside the current directory, will be symlinked into ``data/`` at runtime)
and ``sample_type`` can be one of ``flam-seq``, ``iso-seq``, ``ont-direct``, or
``ont-cdna``. 

Annotation file specifications
------------------------------

Several additional files specify the annotations. 

These include the primary genome assembly (FASTA) and gene annotations (gtf), 
which can be downloaded from e.g. ENSEMBL. 

For Drosophila, the command ``make drosophila_annotation`` will automatically
fetch and download the files used in running the pipeline for the publication 
(Ensembl dm6, build 79). For other organisms, the primary
genome assembly and annotations can be fetched from ENSEMBL, and 
using the 
`UCSC toolkit <https://genome.ucsc.edu/goldenPath/help/twoBit.html>`_ as 
documented will generate the required ``genome.chrom.sizes`` file: 

:: 

   faToTwoBit genome.fa genome.2bit
   twoBitInfo genome.2bit stdout | sort -k2rn > genome.chrom.sizes

These tools are installed into the `cia-snakemake` environment upon
environment creation. 

Custom annotations
------------------

In addition to the canonical annotation files, 
several additional files play a critical role in the correction and 
filtering of long read sequencing assemblies: 

1. ``combined.rds.clusters.new.gff`` contains the highly confident positions 
   of 3´ ends. This file serves as an essential component of the pipeline, 
   as it guides the 3´-UTR correction and filtering process. The generation 
   of the database steps and quality metrics is explained in detailed 
   in this `vignette <https://hilgers-lab.github.io/polyADataBase/docs/polyADatabase.html>`_ . This file can be replaced depending on the organism.  
2. ``promoter.db.edp12.bed`` serves as a reference database indicating the 
   locations of validated transcription start sites, obtained from the EDP Eukaryotic Promoter Database for `Drosophila melanogaster <https://epd.epfl.ch/drosophila/drosophila_database.php?db=drosophila>`_. This file is utilized 
   during the FLAIR assembly step to guide the filtering of isoforms resulting 
   from sequencing artifacts occurring at the 5´ end of the isoforms. 
3. ``splice_junctions_filtered.tab`` Short read sequencing SJ.out files `STAR <https://github.com/alexdobin/STAR>`_. Example file in `here <https://github.com/hilgers-lab/SaiLoR/blob/master/inst/exdata/short_read_junctions.SJ.out.tab>`_. We recommend to pull SJ.out into a single SJ.out from many experiments and filter by min counts. This file is required for FLAIR assembly to 
   correct long reads for possible base mismatches and refine exon boundaries. 
4. ``sqanti.polya.list`` txt file contains a list of poly(A) hexamers that are scanned 
   during the sqanti qc step. This list is relevant for the classification of 
   real 3´ ends of the assemblies.

Annotation locations
--------------------

The default location for the above files is in the ``db`` directory, but they can be 
specified in ``config/config.yaml`` as key: value pairs. For example,

::

   annotation: 
      genome.fasta: /data/repository/organisms/dm6_ensembl/genome_fasta/genome.fa

will symlink from the specified the path to ``annotation/genome.fasta``. 

Pipeline overview
-----------------

#. The first steps are to symlink the data and annotation from those specified 
   during configuration to the ``annotation`` and ``data`` directories 
   respectively. 
#. Minimap and STAR indices are created in the ``index`` directories.
#. Minimap or STARlong are used to map the reads, depending on the
   mapper variable associated with the sample type in ``config/config.yaml``. 
   By default, FLAM-seq and Iso-seq use STARlong (as used in 
   `the FLAM-seq protocol <https://www.nature.com/articles/s41592-019-0503-y>`_)
   whereas the Nanopore protocols use minimap with the `-ax splice -uf` options.
#. The resulting BAM file is indexed. 
#. The BAM file is converted to BED12 format. 
#. Environments for FLAIR and SQANTI are set up manually. Note that SQANTI uses
   the ``sqanti_git_commit`` specified in ``config/config.yaml`` to ensure
   that the pipeline matches the . 
#. FLAIR is used to correct misaligned splice sites from the BED12 file 
   using genome annotations. 
#. FLAIR is used to collapse corrected reads into high-confidence isoforms. 
#. The ``R/RDSFLAM.endGuidedCorrection.R`` is used to correct the ends using
   Nanopore direct RNA-seq or FLAM-seq data. 

Please refer to the STAR Methods in Alfonso-Gonzalez et al. for further details.

End correction
--------------
The end correction script at the end uses the polyadenylation site (PAS) 
database based on FLAM-seq and Nanopore direct RNA-seq data to 
correct the SQANTI-annotated isoforms. The polyadenylation site (PAS)
database reference file is available at ``db/combined.rds.clusters.new.gff`` 
and can be replaced according to the organism or updated with new data. 
The generation of the database steps and quality metrics is explained in 
detailed in this `vignette <https://hilgers-lab.github.io/polyADataBase/docs/polyADatabase.html>`_. 

The end correction script refines and filter individual assemblies 
generated by each method using the polyadenylation site (PAS) database. 
The subsequent filtering parameters were implemented: 

#. All isoforms that overlapped an annotated 3ʹ end within a 100 nt window are retained. 
#. 3ʹ ends identified in the assembly that were more distal than those in the database
   were only retained if they were located within the reference annotation and harbored an 
   AATAAA signal. To correct 3ʹ ends, we initially created 3ʹ untranslated region (UTR) bins 
   using the PAS database, which were established from the end of the open reading frame 
   to the most distal PAS, between each consecutive PAS. 
#. Isoform 3ʹ ends located within the final bin of the 3ʹ UTR (between the two distal-most PASs) 
   were corrected to the most distal bin, provided that the isoform covered greater than 10% of the last bin.

.. 
   _image:: _static/yourimage.png

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
