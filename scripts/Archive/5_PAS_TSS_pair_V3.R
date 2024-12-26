## update: 11-13-2024 
## To remove genes from a GTF file where the distance between two or more of its non-overlapping transcripts exceeds 10 kb
## and remove genes that have transcripts in different strands, such as Gm2004 in mouse refgene.

library(optparse)
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input BAM file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="TSS-exon_coordination.txt",
              help="Output TSS-exon_coordination file", metavar="file"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="Reference annotation gtf file", metavar="file"),
  make_option(c("-p", "--pas"), type="character", default=NULL,
              help="Reference polyA.bed", metavar="file")
)
#Examples: Rscript 5_PAS_TSS_pair_V2.R -i Encode.TSS.TES2.bam -g /Users/yangyalan/results/Long-read-APA-pipeline/reference/hg38.refGene.gtf -p /Users/yangyalan/results/Long-read-APA-pipeline/Encode/Encode.brain.new.PAS.bed -o TSS-TES.coordination.chiqtest.txt
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
bam <- opt$input
pas.bed <- opt$pas
annot_path <- opt$gtf
output <- opt$output


library(dplyr)
library(data.table)
library(tidyr)
library(GenomicRanges)
library(edgeR)
library(stringr)
prepareTSSDatabase <- function(annotation, tss.window) {
  strandSort <- function(x) {
    c(
      GenomicRanges::sort(x[x@strand == "+"], decreasing = FALSE),
      GenomicRanges::sort(x[x@strand == "-"], decreasing = TRUE)
    )
  }
  # Build 5'-3' links data base
  # Exon ids
  txdb <- txdbmaker::makeTxDbFromGRanges(annotation)
  ebt <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  t2g <- AnnotationDbi::select(txdb,
                               keys = names(ebt),
                               keytype = "TXNAME",
                               columns = "GENEID"
  )
  e2 <- BiocGenerics::unlist(ebt)
  e2$transcript_id <- names(e2)
  e2$gene_id <- t2g$GENEID[match(e2$transcript_id, t2g$TXNAME)]
  e2$exon_id <- e2$exon_name
  e2$exon_name <- NULL
  e2$type <- "exon"
  names(e2) <- NULL
  mcols(e2) <- mcols(e2)[, c(
    "exon_id", "exon_rank",
    "transcript_id", "gene_id", "type"
  )]
  bins <- list()
  # TSS data base
  # take first position per transcript and make it single nt
  tss.bins <-
    strandSort(
      plyranges::mutate(
        plyranges::anchor_5p(
          dplyr::filter(e2, exon_rank == 1)), width = 1))
  # make unique TSS starts merging in a 50nt window.
  message("Preparing TSSBase")
  tss.base <-
    strandSort(
      GenomicRanges::makeGRangesFromDataFrame(
        reshape::melt(GenomicRanges::reduce(
          GenomicRanges::split(tss.bins, ~gene_id),
          min.gapwidth = tss.window
        )),
        keep.extra.columns = TRUE
      )
    )
  tss.base <-
    tibble::as_tibble(tss.base) %>%
    dplyr::group_by(value.group_name) %>%
    dplyr::mutate(
      count = paste0(value.group_name, ":P",
                     sprintf("%02d", sequence(dplyr::n())))) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  
  message("TSS database sucessfully created")
  return(tss.base)
}


#annot_path <- "/Users/yangyalan/results/Long-read-APA-pipeline/reference/hg38.refGene.gtf"
message("preparing TSS reference")
refAnnotation <- rtracklayer::import.gff(annot_path)
valid_chromosomes <- paste0("chr", c(1:22, "X","Y"))
refAnnotation <- refAnnotation[seqnames(refAnnotation) %in% valid_chromosomes]
transcripts <- refAnnotation[refAnnotation$type == "transcript"]
gene_strands <- transcripts %>% 
  as.data.frame() %>%
  group_by(gene_id) %>%
  summarize(unique_strands = n_distinct(strand))

# Identify genes that appear on multiple strands
genes_on_multiple_strands <- gene_strands$gene_id[gene_strands$unique_strands > 1]

# Remove these genes from further analysis
transcripts <- transcripts[!transcripts$gene_id %in% genes_on_multiple_strands]

transcript_list <- split(transcripts, transcripts$gene_id)
# Identify genes where all distances between consecutive transcripts are <= 100 kb
genes_to_keep <- sapply(transcript_list, function(gene_transcripts) {
  gene_transcripts <- sort(gene_transcripts)
  distances <- start(gene_transcripts)[-1] - end(gene_transcripts)[-length(gene_transcripts)]
  all(distances <= 10000)
})

valid_genes <- names(genes_to_keep[genes_to_keep])
refAnnotation <- refAnnotation[refAnnotation$gene_id %in% valid_genes]

refExons <- refAnnotation[refAnnotation$type == "exon"]
TSS_database <- prepareTSSDatabase(refExons, tss.window=50)

TSS.anno <- as.data.frame(TSS_database)
write.table(TSS.anno,"TSS.anno.txt",sep="\t",quote=F, row.names = F)

TSS.anno$Tss_id <-  paste0(TSS.anno$seqnames, ":", TSS.anno$start, ":", TSS.anno$strand)
TSS.anno <- TSS.anno %>% dplyr::select(count, Tss_id)
colnames(TSS.anno) <- c("promoter_id","tss_id")

## downstream and upstreamm of the TSSs
#tss_coordinate_base2 <- resize(tss_coordinate_base, width = width(tss_coordinate_base) + 100, fix = "center")
#export(tss_coordinate_base2, "hg38.refGene.TSS.bed", format = "BED")

message("preparing poly(A) sites")
pas <- fread(pas.bed, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pas <- pas[pas$Feature != "intergenic"]
pas <- pas %>% dplyr::filter(gene_count >= 10 & PAU >=0.05 & PAU <= 0.95) %>%  dplyr::select(Chr,Start,End,PAS_ID, gene_name,Strand) %>% arrange(Chr, Start, End)
pas$Start=pas$Start+1
pas.base <- tibble::as_tibble(pas) %>% dplyr::group_by(Chr,Strand,gene_name) %>% 
  dplyr::mutate(count = paste0(gene_name, ":T", sprintf("%02d", sequence(dplyr::n())))) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

pas.anno <- DataFrame(count=elementMetadata(pas.base)$count,
                      PAS_ID_new=elementMetadata(pas.base)$PAS_ID)
write.table(pas.anno,"pas.anno.txt",sep="\t",quote=F, row.names = F)

colnames(pas.anno) <- c("tes_id","pas_id")
pas.anno <- as.data.frame(pas.anno)

message("poly(A) database sucessfully created")

message("reading bam file")
#bam <- "/Users/yangyalan/results/Long-read-APA-pipeline/Encode/TSS-TES/Encode.TSS.TES2.bam"
bamAlignments <- GenomicAlignments::readGAlignments(bam, use.names = TRUE)
alignments <- GenomicRanges::GRanges(bamAlignments)
alignments$name <- names(bamAlignments)
names(alignments) <- NULL


prepareForCountStarts <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  GenomicRanges::end(pos) <- GenomicRanges::start(pos) + window
  GenomicRanges::start(neg) <- GenomicRanges::end(neg) - window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}

prepareForCountEnds <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  GenomicRanges::start(pos) <- GenomicRanges::end(pos) - window
  GenomicRanges::end(neg) <- GenomicRanges::start(neg) + window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}

readTSSassignment <- function(startsAlignemnts, TSSCoordinate.base) {
  tssDb <- TSSCoordinate.base
  ovlps <- GenomicRanges::findOverlaps(startsAlignemnts , tssDb , maxgap = 50)
  promoterIds <- tssDb[subjectHits(ovlps),]$count
  promoterStarts <- GenomicRanges::start(tssDb[subjectHits(ovlps),])
  startsAlignemnts2 <- startsAlignemnts[queryHits(ovlps),]
  startsAlignemnts2$promoter_id <- promoterIds
  startsAlignemnts2$promoterStarts <- promoterStarts
  startsAlignemnts2$dist2Assignment <-
    abs(GenomicRanges::start(startsAlignemnts2) - startsAlignemnts2$promoterStarts)
  startsAlignemnts2 <-
    startsAlignemnts2[order(startsAlignemnts2$dist2Assignment, decreasing = FALSE),]
  startsAlignemnts2 <-
    startsAlignemnts2[!duplicated(startsAlignemnts2$name),]
  return(startsAlignemnts2)
}

readTESassignment <- function(endsAlignements, TESCoordinate.base) {
  tesDb <- TESCoordinate.base
  ovlps <- findOverlaps(endsAlignements , tesDb , maxgap = 150)
  tesIds <- tesDb[subjectHits(ovlps),]$count
  endStarts <- GenomicRanges::start(tesDb[subjectHits(ovlps),])
  endsAlignemnts2 <- endsAlignements[queryHits(ovlps),]
  endsAlignemnts2$tes_id <- tesIds
  endsAlignemnts2$endStarts <- endStarts
  endsAlignemnts2$dist2Assignment <-
    abs(GenomicRanges::start(endsAlignemnts2) - endsAlignemnts2$endStarts)
  endsAlignemnts2 <-
    endsAlignemnts2[order(endsAlignemnts2$dist2Assignment, decreasing = FALSE),]
  endsAlignemnts2 <-
    endsAlignemnts2[!duplicated(endsAlignemnts2$name),]
  return(endsAlignemnts2)
}
countLinks <- function(alignmentsFile, TSSDatabase, PAS) {
  # make single nt starts
  startsAlignemnts <-  prepareForCountStarts(alignmentsFile, 1)
  startAlignments <- readTSSassignment(startsAlignemnts, TSSDatabase)
  startAlignments <-
    startAlignments %>%
    as.data.frame(.) %>%
    dplyr::select(name, promoter_id)
  # make single nt ends
  endsAlignemnts <-  prepareForCountEnds(alignmentsFile, 1)
  endsAlignemnts <- readTESassignment(endsAlignemnts, PAS)
  endsAlignemnts <- endsAlignemnts %>% as.data.frame(.) %>% dplyr::select(name, tes_id)
  # make pairs
  pairsTested <-
    left_join(endsAlignemnts,
              startAlignments,
              by = "name") %>% dplyr::filter(gsub("\\:.*", "", tes_id) == gsub("\\:.*", "", promoter_id)) 
  
  pairsTested <- pairsTested %>% mutate(gene_id =  gsub(":.*", "", pairsTested$promoter_id))

  pairsTested <- pairsTested %>%
    mutate(pairs_id = paste0(
      gene_id,
      ":",
      gsub(".*:", "", pairsTested$promoter_id),
      ":",
      gsub(".*:", "", pairsTested$tes_id)
    ))
  # summarize only reads that expand both pairs
  countedPairs <- pairsTested %>%
    group_by(pairs_id) %>%
    tally() %>%
    dplyr::rename(read_counts = n)
  countedPairsFinal <-
    left_join(
      countedPairs,
      pairsTested %>%
        dplyr::distinct(pairs_id, .keep_all = TRUE) %>%
        dplyr::select(gene_id, pairs_id),
      by = "pairs_id"
    )
  # normalize expression in CPMs
  counts_final <- list()
  dd <- edgeR::DGEList(counts = countedPairsFinal$read_counts)
  dge <- edgeR::calcNormFactors(dd)
  countedPairsFinal$cpm <- as.numeric(edgeR::cpm(dge))
  counts_final$countedPairsFinal <- countedPairsFinal
  counts_final$pairsTested <- pairsTested
  
  return(counts_final)
}

message("counting links between TSS and TES")
countData <- countLinks(alignments, TSS_database, pas.base)

counts <- data.frame(countData$countedPairsFinal)
counts <- counts %>% group_by(gene_id) %>%  filter(dplyr::n() >= 4) %>%  ungroup()
counts <- counts %>% tidyr::separate(pairs_id, into = c("gene_id", "Promoter", "PAS"), sep = ":")
write.table(counts,"TSS-polyA.coordination.count.txt",sep="\t",quote=F, row.names = F)

counts_anno <- counts 
counts_anno$promoter_id <- paste0(counts_anno$gene_id, ":", counts_anno$Promoter)
counts_anno$tes_id <- paste0(counts_anno$gene_id, ":", counts_anno$PAS)
counts_anno <- left_join(counts_anno,TSS.anno, by="promoter_id")
counts_anno <- left_join(counts_anno,pas.anno, by="tes_id")

merged_data1 <- counts_anno %>% select(gene_id,tss_id) %>% unique() %>% group_by(gene_id) %>%
  summarize(tss_id = str_c(tss_id, collapse = ";"))

merged_data2 <- counts_anno %>% select(gene_id,pas_id) %>% unique() %>% group_by(gene_id) %>%
  summarize(pas_id = str_c(pas_id, collapse = ";"))

qname <- countData$pairsTested$name
writeLines(unlist(qname), "TSS-TES.full-length.txt")

gene_matrices <- counts %>%
  group_by(gene_id) %>%
  group_split() %>%
  lapply(function(data) {
    data %>%
      select(Promoter, PAS, read_counts) %>%
      pivot_wider(names_from = Promoter, values_from = read_counts) %>%
      replace(is.na(.), 0)
  })
names(gene_matrices) <- counts$gene_id %>% unique()

message("Peforming statistics analysis")
chisq_results <- lapply(gene_matrices, function(matrix) {
  matrix <- as.matrix(matrix[,-1])
  if (nrow(matrix) >= 2 && ncol(matrix) >= 2) {
    chisq.test(matrix)$p.value
  } else {
    NA  
  }
})
chisq_df <- data.frame(
  gene_id = names(chisq_results),
  p_value = unlist(chisq_results),
  stringsAsFactors = FALSE
)
chisq_df <- chisq_df[order(chisq_df$p_value), ]
chisq_df$FDR <- p.adjust(chisq_df$p_value, method = "BH")
chisq_df$sig <- ifelse(chisq_df$FDR < 0.05, "TRUE", "FALSE")

chisq_df <- left_join(chisq_df,merged_data1, by="gene_id")
chisq_df <- left_join(chisq_df,merged_data2, by="gene_id")

write.table(chisq_df, output, sep="\t",quote=F,row.names = F)
