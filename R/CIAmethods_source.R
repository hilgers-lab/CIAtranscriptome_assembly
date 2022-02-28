#############################################################
## Assembly correction with FLAM-seq/RDS 3'ends source 
## Nov 2020 
#############################################################


suppressMessages(require(rtracklayer))
suppressMessages(require(dplyr))
suppressMessages(require(GenomicFeatures))
suppressMessages(require(ggplot2))
suppressMessages(require(ggpubr))
suppressMessages(require(reshape2))
suppressMessages(require(GenomicFeatures))
suppressMessages(require(GenomicRanges))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(tximport))
suppressMessages(require(Biostrings))
suppressMessages(require(BSgenome.Dmelanogaster.UCSC.dm6))
suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))
suppressMessages(require(ggpubr))
suppressMessages(require(RColorBrewer))
suppressMessages(require(data.table))
suppressMessages(require(UpSetR))
suppressMessages(require(tidyverse))


##configs 
hexamers <- list()   #output list of examers

#APA motifs  
APAmotifs <- c("AATAAA", "ATTAAA","AATATA", "AAGAAA",
               "AATACA", "AATAGA", "AATATA", "AATGAA", 
               "ACTAAA", "CATAAA", "GATAAA", "TATAAA", "TTTAAA" ) 

## Functions -------------------------------------------------------------------------------------
strandAwareFindOverlapEnd <- function(queryData, annotation, maxgap) {
  hitpos <- findOverlaps(queryData[queryData@strand=="+"], 
                         annotation[annotation@strand=="+"], 
                         type="end",
                         maxgap = maxgap)
  hitneg <- findOverlaps(queryData[queryData@strand=="-"], 
                         annotation[annotation@strand=="-"], 
                         type="start",
                         maxgap = maxgap)
  
  bin.exon <- c( queryData[queryData@strand=="+"][unique(queryHits(hitpos))] , 
                 queryData[queryData@strand=="-"][unique(queryHits(hitneg))] ) 
  
  return(bin.exon)
}

getLastExon <- function(seq.gr){
  seq.gr.pos <- sort.GenomicRanges(seq.gr[seq.gr@strand=="+"], by= ~end, decreasing = T) #the first is the last 
  seq.gr.neg <- sort.GenomicRanges(seq.gr[seq.gr@strand=="-"], by= ~start, decreasing = F) #the first is the last 
  seq.gr <- c(seq.gr.pos, seq.gr.neg)
  seq.gr.ends <- seq.gr[!duplicated(seq.gr$gene_id),]
  return(seq.gr.ends)
} 

getLastExonperTranscript <- function(seq.gr){
  seq.gr.pos <- sort.GenomicRanges(seq.gr[seq.gr@strand=="+"], decreasing = T) #the first is the last 
  seq.gr.neg <- sort.GenomicRanges(seq.gr[seq.gr@strand=="-"], decreasing = F) #the first is the last 
  seq.gr <- c(seq.gr.pos, seq.gr.neg)
  seq.gr.ends <- seq.gr[!duplicated(seq.gr$transcript_id),]
  return(seq.gr.ends)
}

getTranscriptsWithDistalExon <- function(seq.gr){
  seq.gr.pos <- sort.GenomicRanges(seq.gr[seq.gr@strand=="+"], by= ~end, decreasing = T) #the first is the last 
  seq.gr.neg <- sort.GenomicRanges(seq.gr[seq.gr@strand=="-"], by= ~start, decreasing = F) #the first is the last 
  seq.gr <- c(seq.gr.pos, seq.gr.neg)
  seq.gr.ends <- seq.gr[!duplicated(seq.gr$gene_id),]
  seq.gr.transcriptswlastexon<- getLastExonperTranscript(seq.gr)
  overlapping.ends <- strandAwareFindOverlapEnd(seq.gr.transcriptswlastexon, seq.gr.ends, maxgap = 3000)
  transcriptsWithOvends <- seq.gr[seq.gr$transcript_id %in% overlapping.ends$transcript_id] 
  transcriptswlastexon<- getLastExonperTranscript(transcriptsWithOvends)
  return(transcriptswlastexon)
} 

createScaffoldBlocks <- function(isoforms.sqanti, LongEnds){ 
  ## function to create scaffold blocks using sqanti ORF prediction and High quality 3'neds from FLAM-seq/RDS 
  ## get stop codons from sqanti class file 
  # pos genes 
  pos.strand  <- isoforms.sqanti %>%
    dplyr::select("associated_gene","CDS_genomic_end", "CDS_genomic_start", "strand", "associated_transcript") %>% 
    filter(strand=="+") %>% 
    group_by(associated_gene) %>% 
    arrange(desc(CDS_genomic_end)) %>% 
    dplyr::slice(1) %>% 
    arrange(desc(CDS_genomic_end))
  
  # negative  strand genes 
  
  neg.strand  <- isoforms.sqanti %>%
    dplyr::select("associated_gene","CDS_genomic_end", "CDS_genomic_start", "strand", "associated_transcript") %>% 
    filter(strand=="-") %>% 
    group_by(associated_gene) %>% 
    arrange(CDS_genomic_end) %>% 
    dplyr::slice(1) %>% 
    arrange(CDS_genomic_end)
  
  ## assign stop codon to 3'ends 
  
  LongEnds.short <- LongEnds
  
  mneg <- makeGRangesFromDataFrame( 
    merge(LongEnds.short[LongEnds.short$strand=="-",], 
          neg.strand, 
          by.x="gene_id", 
          by.y="associated_gene", 
          all.x=TRUE),  keep.extra.columns = TRUE)
  strand(mneg) <- mneg$strand.x
  
  message("Build from sqanti stop codon database with 3'ends")
  
  
  mneg$dif  <- start(mneg) - mneg$CDS_genomic_end  
  
  message("Removing: ",length(mneg[-which(mneg$dif<0),]), " transcript 3'ends that are before stop codon in negative strand ")
  
  #genes.with.apa.alternative.stop.splicing. <- mneg[-which(mneg$dif<0),]
  mneg <- mneg[which(mneg$dif<0),]
  
  end(mneg) <-  mneg$CDS_genomic_end
  
  # strand(mneg) <- mneg$strand
  
  ### Disjoin 3'UTR to generate 3'UTR block bins in a gene specific manner  
  
  ## creation of 3'UTR scaffolds blocks 
  
  ## Now for positive strand. 
  
  mpos <- makeGRangesFromDataFrame( 
    merge(LongEnds.short[LongEnds.short$strand=="+",], 
          pos.strand, 
          by.x="gene_id", 
          by.y="associated_gene", 
          all.x=TRUE), keep.extra.columns = TRUE)
  
  strand(mpos) <- mpos$strand.x
  mpos$dif  <- start(mpos) - mpos$CDS_genomic_end  
  #genes.with.apa.alternative.stop.splicing. <- mpos[which(mpos$dif<0),]
  message("Removing: ",length(mpos[which(mpos$dif<0),]), " transcript 3'ends that are before stop codon in negative positive ")
  mpos <- mpos[which(mpos$dif>0),]
  
  message("extend 3'end to stop codon coordinate")
  start(mpos) <-  mpos$CDS_genomic_end
  ##  Disjoin extended 3'ends in 3'UTR scaffold blocks. 
  ## creation of 3'UTR scaffolds blocks  
  message("disjoin 3'ends and make scaffold blocks")
  mneg.gr.blocks <- unlist(disjoin(split( mneg, ~gene_id)))
  mnpos.gr.blocks <- unlist(disjoin(split( mpos, ~gene_id)))
  
  blocks <- c(mneg.gr.blocks, mnpos.gr.blocks)
  
  return(blocks)
  
}



blockGuidedcorrection <- function(blocks, isoform.assembly, LongWindow, shortWindow ){
  ## create object that contains intermediate steps 
  corrected.res <- list()
  
  ## add gene id 
  blocks$gene_id <- names(blocks)
  names(blocks) <- NULL
  
  ## compare  distal to 3'ends from FLAM-seq/RDS 
  distalExon <- getLastExon(blocks) # comparison only to most distal block 
  ## get 3'UTRs per gene 
  distalExonTx <- getTranscriptsWithDistalExon(isoform.assembly) #comparison to all 3'UTRs  ##i made a big change here replaced getTranscriptsLastExonw with getTranscriptsWithDistalExon
  
  diexon.ov <- findOverlaps(distalExonTx, 
                            distalExon)
  
  
  df.data <- cbind( as.data.frame(distalExonTx[queryHits(diexon.ov)]) ,  
                    as.data.frame(distalExon[subjectHits(diexon.ov)])  )  
  
  
  print( nrow(df.data) ) 
  colnames(df.data) <- make.unique(colnames(df.data))
  
  message("Identify 3'ends for correction")
  
  df.data$dis2Cleavage <- ifelse(df.data$strand=="+", df.data$end - df.data$end.1, 
                                 df.data$start - df.data$start.1)
  
  df.data$dis2start <- ifelse(df.data$strand=="+", df.data$end - df.data$start.1, 
                              df.data$start - df.data$end.1)
  
  
  df.data$end_quality <- ifelse(df.data$strand=="+" & df.data$dis2Cleavage<c(-shortWindow) |
                                  df.data$strand=="-" & df.data$dis2Cleavage>shortWindow, "shorter", "other")
  
  df.data$end_quality <- ifelse(df.data$strand=="+" & df.data$dis2Cleavage>LongWindow|
                                  df.data$strand=="-" & df.data$dis2Cleavage<c(-LongWindow), "longer"  , df.data$end_quality )
  
  #ifelse(res.correction$ends.data$strand=="+", res.correction$ends.data$dis2Cleavage /  res.correction$ends.data$width.1 , res.correction$ends.data$dis2start /   res.correction$ends.data$width.1) 
  
  # genes with single 3'ends are not corrected to avoid over extension
  
  
  single.end.gene.id <- as.data.frame(blocks) %>% 
    group_by(gene_id) %>% 
    tally() %>% 
    filter(n<2) %>% ## warning here
    pull(1)
  
  df.data$ends_per_gene <- ifelse(df.data$gene_id %in% single.end.gene.id, "single.end", "more")
  
  corrected.res$ends.data <-  df.data 
  # transcripts can be corrected if:
  ## they have more than 1 FLAM-seq end reported 
  ## the assembly is shorter than the distal 3'end reference but theres not a longer one
  print( nrow(df.data) )
  
  transcripts2Correct <- df.data %>% 
    filter(!ends_per_gene == "single.end") %>%
    group_by(gene_id) %>%
    filter(!any(end_quality== "longer")) %>%
    ungroup() %>%
    filter(end_quality== "shorter")
  
  #new: correct only those genes in which the shortest covers at least 25% of the most distal scaffold block. 
  transcripts2Correct$prop.coverage <- ifelse(transcripts2Correct$strand=="-", 
                                                   abs(transcripts2Correct$dis2start) / transcripts2Correct$width.1 , 
                                                   abs(transcripts2Correct$dis2start) / transcripts2Correct$width.1   ) 
  
  transcripts2Correct <- subset(transcripts2Correct, prop.coverage>0.1)
  
  
  message("LongEnd read correction for: ",  nrow(transcripts2Correct), " transcripts")
  
  
  #create unique 3'UTR distal identifier 
  
  transcripts2Correct$neu.id  <- paste(transcripts2Correct$transcript_id, transcripts2Correct$end, sep=".")
  
  transcripts2Correct <- transcripts2Correct %>% dplyr::select("neu.id", "start.1", "end.1", "end_quality", "gene_id.1")
  
  corrected.res$corrected.ends <- transcripts2Correct 
  
  
  ## create unique id for all exons of reference annotation
  
  assembly.df <- as.data.frame(isoform.assembly)
  
  assembly.df$neu.id  <- paste(assembly.df$transcript_id, assembly.df$end, sep=".")
  
  assembly.df <- merge(assembly.df, 
                       transcripts2Correct, 
                       by="neu.id", 
                       all.x=TRUE)
  
  assembly.df$end_quality <- ifelse(is.na(assembly.df$end_quality) , 
                                    "other.ends",  assembly.df$end_quality )
  
  
  ## strand aware correction only for those in which the 3'end is shorter than FLAM-seq and that falls in the block bin those 3'ends are corrected. 
  pos <- assembly.df[assembly.df$strand=="+",]
  neg <- assembly.df[assembly.df$strand=="-",]
  
  pos$end <- ifelse(pos$end_quality=="shorter", 
                    pos$end.1, pos$end )
  neg$start <-ifelse(neg$end_quality=="shorter", 
                     neg$start.1, neg$start )
  
  assembly.df.neu <- rbind(pos, neg) 
  
  message("Sucessfully corrected transcript 3'ends: ",  nrow(transcripts2Correct), "  transcripts") 
  assembly.df.neu <- assembly.df.neu[2:13]
  corrected.res$extended.assembly <- makeGRangesFromDataFrame(assembly.df.neu, keep.extra.columns = T)
  
  return(corrected.res)
  
  
}

classifyLonger <- function(ends.data, ens.annot){ 
  df.long <- subset(ends.data, end_quality=="longer") %>% dplyr::select("seqnames", "start", "end", "strand", 
                                                                        "source", "gene_id", "transcript_id", "dis2Cleavage", 
                                                                        "end_quality", "ends_per_gene")
  df.long.gr <- makeGRangesFromDataFrame(df.long, keep.extra.columns = T)
  ## get most distal from the longer than reference 
  most.distal <- getTranscriptsWithDistalExon(df.long.gr)
  df.long.gr$utr_class <- ifelse(df.long.gr$transcript_id %in% most.distal$transcript_id, "ultra.distal", "other")
  most.distal2 <- getTranscriptsWithDistalExon(df.long.gr[-which(df.long.gr$utr_class=="ultra.distal"),])
  df.long.gr$utr_class <- ifelse(df.long.gr$transcript_id %in% most.distal2$transcript_id, "second.distal", df.long.gr$utr_class)
  df.long.gr.kept.gr.utr <- gRanges2Biostrings3UTR(df.long.gr)
  dfs <- annotatepolyAsignals(df.long.gr.kept.gr.utr, APAmotifs, hexamers)
  df.long.gr.signals <- merge(as.data.frame(df.long.gr.kept.gr.utr)[c("transcript_id", "utr_class", "dis2Cleavage", "end_quality", "ends_per_gene")],
                              dfs, by="transcript_id", all.x=TRUE)
  
  ### check if it match ensembl annotation 
  utr.bin <- ens.annot[ens.annot$type=="three_prime_utr"]
  utr.bin.short <- gRanges2Biostrings3UTR(utr.bin)
  df.long.gr.short <- gRanges2Biostrings3UTR(df.long.gr)
  stringent.in.ens <- as.data.frame( strandAwareFindOverlapEnd(df.long.gr.short , utr.bin.short, 100) ) %>% 
    dplyr::select("transcript_id") %>% pull() 
  df.long.gr.signals$is.Ens.End <- ifelse(df.long.gr.signals$transcript_id %in% stringent.in.ens, TRUE,FALSE)
  return(df.long.gr.signals)
  
  
}

gRanges2Biostrings3UTR <- function(grset){
  pos <- grset[grset@strand=="+"]
  start(pos) <- end(pos)-50
  neg <- grset[grset@strand=="-"]
  end(neg) <- start(neg)+50
  grset <- c(pos,neg)
  grset <- as.data.frame(grset)
  grset$seqnames <- paste0("chr", grset$seqnames)
  x <- subset(grset, seqnames %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY") )
  link_on.gr <- makeGRangesFromDataFrame(x, keep.extra.columns = T)
  return(link_on.gr)
}

# Signals and cluster stats 
gRanges2Biostrings3UTR <- function(grset){
  pos <- grset[grset@strand=="+"]
  start(pos) <- end(pos)-50
  neg <- grset[grset@strand=="-"]
  end(neg) <- start(neg)+50
  grset <- c(pos,neg)
  grset <- as.data.frame(grset)
  grset$seqnames <- paste0("chr", grset$seqnames)
  x <- subset(grset, seqnames %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY") )
  link_on.gr <- makeGRangesFromDataFrame(x, keep.extra.columns = T)
  return(link_on.gr)
}

getLastExonperTranscript <- function(seq.gr){
  seq.gr.pos <- sort.GenomicRanges(seq.gr[seq.gr@strand=="+"], decreasing = T) #the first is the last 
  seq.gr.neg <- sort.GenomicRanges(seq.gr[seq.gr@strand=="-"], decreasing = F) #the first is the last 
  seq.gr <- c(seq.gr.pos, seq.gr.neg)
  seq.gr.ends <- seq.gr[!duplicated(seq.gr$transcript_id),]
  seq.gr.ends
} 

getPolyAsignal <- function(APAmotifs, hexamers, seqs){
  if(rlang::is_empty(APAmotifs)){ # no motifs remaining
    return(hexamers)
  } else {
    apaMotif = APAmotifs[1] #get first entry in list
    dnaseq <- DNAString(apaMotif) #canonical 1
    hexamers[[apaMotif]] <-unlist(vmatchPattern(dnaseq, seqs))
    seqs <- seqs[!names(seqs)  %in% names( hexamers[[apaMotif]] ),]
    APAmotifs <- APAmotifs[!APAmotifs %in%  apaMotif] # remove motif from list
    getPolyAsignal(APAmotifs, hexamers, seqs)
  }
}

getDistance2poly<- function(dataframe){
  link_onp.gr <- makeGrangesFromThreeLinks(dataframe)
  #make short fragments with 50'nt windows 
  link_on.gr.short <- gRanges2Biostrings3UTR(link_onp.gr) #outputs granges
  #make seqs
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, link_on.gr.short)
  names(seqs) <- link_on.gr.short$name
  signals <- getPolyAsignal(APAmotifs, hexamers,seqs)
  signals2 <- lapply(signals, as.data.frame)
  signals2 <- bind_rows(signals2, .id = "column_label")
  signals2$dist_polya  <- signals2$end-50
  canonical <- c("AATAAA", "ATTAAA","AATATA")
  signals2$signal_short <- ifelse(signals2$column_label %in% canonical, 
                                  signals2$column_label, "non_canonical" )
  return(signals2)
}

getSignalDist <- function(x) {
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, x)
  names(seqs) <- x$name
  signals <- getPolyAsignal(APAmotifs, hexamers,seqs)
  signals2 <- lapply(signals, as.data.frame)
  signals2 <- bind_rows(signals2, .id = "column_label")
  signals2$dist_polya  <- signals2$end-50
  canonical <- c("AATAAA", "ATTAAA","AATATA")
  signals2$signal_short <- ifelse(signals2$column_label %in% canonical, 
                                  signals2$column_label, "non_canonical" )
  return(signals2)
  
}

relativeNucleotidePolyAsite <- function(grs){
  short_utrs <- gRanges23UTR5nt(grs)
  flank = 100
  pas.flank <- flank(short_utrs , width = flank, both = TRUE)
  pas.flank_seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, pas.flank)
  prob <- consensusMatrix(pas.flank_seq,as.prob = T)
  df <- melt(prob)
  prob1 <- prob 
  prop3 <- list(prob1, prob)
  df <- melt(prop3)
  levels(df$Var1)[4] <- 'U'
  df = subset(df, grepl("A|C|G|U", Var1))
  return(df) 
}


gRanges23UTR5nt <- function(grset){
  pos <- grset[grset@strand=="+"]
  start(pos) <- end(pos)-1
  neg <- grset[grset@strand=="-"]
  end(neg) <- start(neg)+1
  grset <- c(pos,neg)
  grset <- as.data.frame(grset)
  #grset$seqnames <- paste0("chr", grset$seqnames)
  #x <- subset(grset, grepl("chr3|chr2|chr4|X", seqnames))
  x <-grset
  link_on.gr <- makeGRangesFromDataFrame(x, keep.extra.columns = T)
  return(link_on.gr)
}

flank = 99

annotatepolyAsignals <- function(gr.object, APAmotifs, hexamers) { 
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, gr.object)
  names(seqs) <- gr.object$transcript_id
  signals <- getPolyAsignal(APAmotifs, hexamers,seqs)
  signals2 <- lapply(signals, as.data.frame)
  signals2 <- bind_rows(signals2, .id = "column_label")
  signals2$dist_polya  <- signals2$end-50
  canonical <- c("AATAAA", "ATTAAA","AATATA")
  signals2$signal_short <- ifelse(signals2$column_label %in% canonical, 
                                  signals2$column_label, "non_canonical" )
  signals.all <- merge(as.data.frame(gr.object)[c("gene_id", "transcript_id")], 
                       signals2, by.x="transcript_id", by.y="names", all.x=TRUE)
  signals.all$polya.signals <- ifelse(is.na(signals.all$column_label), 
                                      "not_detected", 
                                      signals.all$column_label)
  signals.all$signal_short <- ifelse(signals.all$polya.signals %in% c("not_detected"), 
                                     "not_detected", signals.all$signal_short )
  return(signals.all) 
  
}

hexamers <- list()


getNucleotideDist <- function(dataset) {
  stringent.removed.gr <- makeGRangesFromDataFrame(dataset[2:13], keep.extra.columns = T)
  stringent.removed.gr.utr <- getLastExonperTranscript(stringent.removed.gr)
  stringent.removed.gr.utr <- gRanges2Biostrings3UTR(stringent.removed.gr.utr)
  stringent.removed.gr.signals <- annotatepolyAsignals(stringent.removed.gr.utr, APAmotifs, hexamers)
  removed.signals.nt <- relativeNucleotidePolyAsite(stringent.removed.gr.utr)
  removed.signals.nt
}


gRanges23UTR100nt <- function(grset){
  pos <- grset[grset@strand=="+"]
  start(pos) <- end(pos)-100
  neg <- grset[grset@strand=="-"]
  end(neg) <- start(neg)+100
  grset <- c(pos,neg)
  grset
}

