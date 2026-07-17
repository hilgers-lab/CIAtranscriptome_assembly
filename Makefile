.PHONY: test run env clean-env drosophila_annotation
ZENODO_URL=https://zenodo.org/record/7438383/files/test.tar.gz
DM6_GENOME_URL=ftp://ftp.ensembl.org/pub/release-79/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz
DM6_GTF_URL=ftp://ftp.ensembl.org/pub/release-96/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gtf.gz
CONDA_ENV=snakemake-cia
CONDA_PACKAGES=snakemake ucsc-twobitinfo ucsc-fatotwobit

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

test: test/flam-seq/.makeflag drosophila_annotation .env-created
	./test.sh
	
test/flam-seq/.makeflag: test.tar.gz
	tar -xzvf $<
	touch $@

test.tar.gz:
	wget https://zenodo.org/record/7438383/files/test.tar.gz

env: .env-created

clean-env:
	rm -f .env-created
	conda env remove -n $(CONDA_ENV)

.env-created: 
	conda create -y -n $(CONDA_ENV) -c conda-forge -c bioconda $(CONDA_PACKAGES)
	touch $@

run:
	./run.sh

drosophila_annotation: db/dm6_ensembl/genome.chrom.sizes db/dm6_ensembl/release-96-genes.gtf

db/dm6_ensembl/genome.fa: 
	mkdir -p $(@D)
	wget -O - $(DM6_GENOME_URL) | gunzip -c > $@

db/dm6_ensembl/release-96-genes.gtf: 
	mkdir -p $(@D)
	wget -O - $(DM6_GTF_URL) | gunzip -c > $@

db/dm6_ensembl/genome.chrom.sizes: db/dm6_ensembl/genome.2bit .env-created
	$(CONDA_ACTIVATE) $(CONDA_ENV) && twoBitInfo $< stdout | sort -k2rn > $@

db/dm6_ensembl/genome.2bit: db/dm6_ensembl/genome.fa .env-created
	$(CONDA_ACTIVATE) $(CONDA_ENV) && faToTwoBit $< $@

