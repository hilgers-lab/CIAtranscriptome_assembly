.. CIA Transcriptome Assembly documentation master file, created by
   sphinx-quickstart on Wed Feb 22 16:03:05 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CIA Transcriptome Assembly's documentation!
======================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Software/hardware requirements
------------------------------

The CIA Transcriptome Assembly has been tested on CentOS 7 and Ubuntu 22.04 LTS
and should work on most Linux64 platforms. At least 32GB of memory are required
for the STAR alignment step. 

Download ``conda`` from 
`here <https://docs.conda.io/en/latest/miniconda.html>`_. Then execute
``make env`` from the base directory to make the base snakemake environment. 
This will result in an environment called ``snakemake-cia``. If you need
to remove the conda environment, execute ``make clean-env``. 

All other software will be installed via conda or during the runtime of the 
pipeline. 

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
---------------------
Several additional files specify the annotations. 

These include the primary genome assembly (FASTA) and gene annotations (gtf), 
which can be downloaded from e.g. ENSEMBL. 

For Drosophila, the command ``make drosophila_annotation`` will automatically
fetch and download the files used in running the pipeline for the publication 
(Ensembl dm6, build 79). For other organisms, using the 
`UCSC toolkit <https://genome.ucsc.edu/goldenPath/help/twoBit.html>`_ as 
documented will generate the ``genome.chrom.sizes`` file: 

``faToTwoBit genome.fa genome.2bit``
``twoBitInfo genome.2bit stdout | sort -k2rn > genome.chrom.sizes``

In addition to the canonical annotation files, several additional files are
used . The files corresponding to Drosophila are included in the ``db/`` folder
in the repository, but will need to be created for other organisms/genome builds.  

#. combined.rds.clusters.new.gff
#. promoter.db.edp12.bed
#. splice_junctions_filtered.tab
#. sqanti.polya.list

The default location for these is in the ``db`` directory, but they can be 
specified . 

combined.rds.clusters.new.gff

Pipeline overview
-----------------

#. The first steps are to symlink the data and annotation from those specified .
   during configuration to the ``annotation`` and ``data`` directories 
   respectively. 
#. Minimap and STAR indices are created in the ``index`` directories.
#. Minimap or STARlong are used to map the reads, depending on the
   mapper variable associated with the sample type in ``config/config.yaml``. 
   By default, FLAM-seq and Iso-seq use STARlong (as used in 
   `the FLAM-seq protocol <https://www.nature.com/articles/s41592-019-0503-y>`_)
   whereas the Nanopore protocols use minimap with the `-ax splice -uf` options.
#. The resulting BAM file is indexed. 
#. Environments for FLAIR and SQANTI are set up manually. Note that SQANTI uses
   the ``sqanti_git_commit`` specified in ``config/config.yaml`` to ensure
   that the pipeline matches the . 


Rule DAG
--------

.. image:: _static/dag.png



Please refer to the STAR Methods in Alfonso-Gonzalez et al. for further details.


Pipeline configuration
----------------------

Testing
-------

A test set is available at `Zenodo <https://doi.org/10.5281/zenodo.7438383>`_,
consisting of one FASTQ for each sample type. 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
