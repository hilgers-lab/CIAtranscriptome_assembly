# Correct long read assemblies using 3'end database 
#################################################################
suppressPackageStartupMessages(library("optparse"))
source("CIAmethods_source.R")
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
option_list <- list( 
  make_option(c("-a", "--annot"), type="character", 
              help="annotation file from assembly in gtf"),
  make_option(c("-e", "--endsRef"), type="character",
              help = "3'end reference database"),
  make_option(c("-s", "--isoClass"), type="character",
              help = "sqanti output *_classification.txt with isoform info"),
  make_option(c("-r", "--refannot"), type="character",
              help = "ref annotation in ensembl format"),
  make_option(c("-o", "--outDir"), type="character",
              help = "ref annotation in ensembl format"),
  make_option(c("-p", "--prefix"), type="character",
              help = "ref annotation in ensembl format")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)
options(stringsAsFactors = FALSE)

gRanges23UTR100nt <- function(grset){
  pos <- grset[grset@strand=="+"]
  start(pos) <- end(pos)-100
  neg <- grset[grset@strand=="-"]
  end(neg) <- start(neg)+100
  grset <- c(pos,neg)
  grset
}


annotFile <-  "/data/processing1/alfonso/rerunAssembly/cdna_seq/cia_isoforms/cia_assemblies/heads_cia_assembly/cia_heads.all.seqs.gtf"
LongEndsFile <-"data/cia_isoforms/clusters.rds/combined.rds.clusters.new.gff"
sqantiClass <- "/data/processing1/alfonso/rerunAssembly/cdna_seq/cia_isoforms/cia_assemblies/heads_cia_assembly/sqanti/cia_heads.all.seqs.sqanti_classification.txt"
ensAnnotFile <- "/data/repository/organisms/dm6_ensembl/ensembl/release-96/genes.gtf"

## Load data  ----------------------------------------------------------------------------------------------
LongEnds <- gRanges23UTR100nt(import.gff(opt$endsRef))
LongEnds <- tibble( as.data.frame(sort( LongEnds, decreasing = TRUE)))
isoforms.sqanti <- readr::read_tsv(opt$isoClass)
uncor.annot <- import.gff(opt$annot)
uncor.annot <- uncor.annot[uncor.annot$type=="exon"]
#ens annotation 
ensAnnot <- import.gff(opt$refannot)
## Generate 3'UTR scafold blocks -------------------------------------------------------------------------------------
# Get stop codons 
## Note: Here as a simplification I use the most distal stop codon to the 3'UTR. 
blocks <- createScaffoldBlocks(isoforms.sqanti, LongEnds) 
## Scaffold guided 3'end correction -------------------------------------------------------------------------------------
### keep only most distal block 
res.correction <- blockGuidedcorrection(blocks, uncor.annot, 300, 50)
### Prepare filtering features  ---------------------------------------------
most.distal.ends  <- classifyLonger(res.correction$ends.data, ensAnnot)
# Classify longer than reference transcripts to keep  
# For ultra distal keep those with cannonical AATAAA OR ensembl annotated 
transcriptsUtraLong <- subset(most.distal.ends, utr_class=="ultra.distal" & signal_short %in% c("AATAAA") | 
                                utr_class=="ultra.distal" & is.Ens.End=="TRUE"  ) %>% dplyr::select(transcript_id)
transcriptsSecondDistal<- subset(most.distal.ends, utr_class=="second.distal" & signal_short%in%c("AATAAA")) %>% 
  dplyr::select(transcript_id)
longToSelect <- c(transcriptsUtraLong$transcript_id, transcriptsSecondDistal$transcript_id )
## Filtering ------------------------------------------------------------------------------------- 
hqLongEnds <- import.gff(opt$endsRef)
uncor.annot <- import.gff(opt$annot)
#this is very important because this ensures the comparison of the most distal 3'end of a transcript against the reference 3'ends
uncor.txs <- uncor.annot[uncor.annot$type=="transcript"]
uncor.annot <- uncor.annot[uncor.annot$type=="exon"]
heads.true.ends.strand.aware <- strandAwareFindOverlapEnd(uncor.txs, hqLongEnds , 250) ##check this in IGV
heads.true.ends.strand.aware <- c(heads.true.ends.strand.aware$transcript_id, 
                                  sub('^([^.]+.[^.]+).*', '\\1', res.correction$corrected.ends$neu.id))
#transcripts_with.known.end.strand.aware  <- unique(heads.true.ends.strand.aware)
transcripts_with.known.end.strand.aware  <- heads.true.ends.strand.aware
transcriptsToKeep <- c(longToSelect , transcripts_with.known.end.strand.aware)
#genesOnlyInassembly <- setdiff(heads.df.neu$gene_id, five.read.ends$name) 
corrected.filtered.assembly <- subset(res.correction$extended.assembly, transcript_id %in% transcriptsToKeep ) 
corrected.removed.assembly <- subset(res.correction$extended.assembly, !transcript_id %in% transcriptsToKeep ) 
utr.from.removed.elements <- getLastExonperTranscript(corrected.removed.assembly)
## export 
export.gff2(corrected.filtered.assembly, paste0(opt$outDir, "/", opt$prefix, ".end.corrected.assembly.gtf") ) 



