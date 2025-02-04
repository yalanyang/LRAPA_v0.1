## update: 11-13-2024 
## To remove genes from a GTF file where the distance between two or more of its non-overlapping transcripts exceeds 10 kb
## and remove genes that have transcripts in different strands, such as Gm2004 in mouse refgene.

library(optparse)
library(AnnotationDbi)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input BAM file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="TSS-exon_coordination.txt",
              help="Output TSS-exon_coordination file", metavar="file"),
  make_option(c("-s", "--tss"), type="character", default=NULL,
              help="Reference annotation gtf file", metavar="file"),
  make_option(c("-p", "--pas"), type="character", default=NULL,
              help="Reference polyA.bed", metavar="file")
)
#Examples: Rscript 5_PAS_TSS_pair_V4.R -i encode.mapping.unique.filter.bam -s /Users/yangyalan/results/Long-read-APA-pipeline/reference/Human.encode.V40.TSS.ref.txt -p /Users/yangyalan/results/Long-read-APA-pipeline/Encode/Encode.brain.new.PAS.bed -o TSS-TES.encode.coordination.chiqtest.txt
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
bam <- opt$input
pas.bed <- opt$pas
tss.anno <- opt$tss
output <- opt$output

library(dplyr)
library(data.table)
library(tidyr)
library(GenomicRanges)
library(edgeR)
library(stringr)

message("loading transcriptional start sites")
TSS.anno<- fread(tss.anno, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

TSS_database <- GRanges(
  seqnames = TSS.anno$seqnames,
  ranges = IRanges(start = TSS.anno$start, end = TSS.anno$end),
  strand = TSS.anno$strand,
  gene_name = TSS.anno$gene_name,
  count = TSS.anno$count
)

TSS.anno$Tss_id <-  paste0(TSS.anno$seqnames, ":", TSS.anno$start, ":", TSS.anno$strand)
TSS.anno <- TSS.anno %>% dplyr::select(count, Tss_id)
colnames(TSS.anno) <- c("promoter_id","tss_id")

message("preparing poly(A) sites")
pas <- fread(pas.bed, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pas <- pas[pas$Feature != "intergenic"]
pas <- pas %>% dplyr::filter(gene_count >= 10 & PAU >=0.05 & PAU <= 0.95) %>%  dplyr::select(Chr,Start,End,PAS_ID, gene_name,Strand) %>% arrange(Chr, Start, End)
pas$Start=pas$Start+1
pas.base <- tibble::as_tibble(pas) %>% dplyr::group_by(Chr,Strand,gene_name) %>% 
  dplyr::mutate(count = paste0(gene_name, ":T", sprintf("%02d", sequence(dplyr::n())))) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

pas.anno <- DataFrame(count=elementMetadata(pas.base)$count,
                      PAS_ID=elementMetadata(pas.base)$PAS_ID)
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
  shortsend <- c(pos, neg)
  return(shortsend)
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
  pairsTested <- left_join(endsAlignemnts, startAlignments, by = "name") 

  pairsTested <- pairsTested %>% dplyr::filter(gsub("\\:.*", "", tes_id) == gsub("\\:.*", "", promoter_id)) 
  
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
  #dd <- edgeR::DGEList(counts = countedPairsFinal$read_counts)
  #dge <- edgeR::calcNormFactors(dd)
  #countedPairsFinal$cpm <- as.numeric(edgeR::cpm(dge))
  counts_final$countedPairsFinal <- countedPairsFinal
  counts_final$pairsTested <- pairsTested
  
  return(counts_final)
}

message("counting links between TSS and TES")
countData <- countLinks(alignments, TSS_database, pas.base)

counts <- data.frame(countData$countedPairsFinal)
counts <- counts %>% tidyr::separate(pairs_id, into = c("gene_id", "Promoter", "PAS"), sep = ":") %>% filter(read_counts > 5)

gene_summary <- counts %>%
  group_by(gene_id) %>%
  summarize(
    num_promoters = n_distinct(Promoter),
    num_PAS = n_distinct(PAS)
  )

genes_with_multiple <- gene_summary %>%
  filter(num_promoters > 1, num_PAS > 1) %>%
  pull(gene_id)

counts <- counts %>%
  filter(gene_id %in% genes_with_multiple) 


#counts <- counts %>% group_by(gene_id) %>%  filter(dplyr::n() >= 4) %>%  ungroup()
write.table(counts,"TSS-polyA.coordination.count.txt",sep="\t",quote=F, row.names = F)

counts_anno <- counts 
counts_anno$promoter_id <- paste0(counts_anno$gene_id, ":", counts_anno$Promoter)
counts_anno$tes_id <- paste0(counts_anno$gene_id, ":", counts_anno$PAS)
counts_anno <- left_join(counts_anno,TSS.anno, by="promoter_id")
counts_anno <- left_join(counts_anno,pas.anno, by="tes_id")

merged_data1 <- counts_anno %>% dplyr::select(gene_id,tss_id) %>% unique() %>% group_by(gene_id) %>%
  summarize(tss_id = str_c(tss_id, collapse = ";"))

merged_data2 <- counts_anno %>% dplyr::select(gene_id,pas_id) %>% unique() %>% group_by(gene_id) %>%
  summarize(pas_id = str_c(pas_id, collapse = ";"))

qname <- countData$pairsTested$name
writeLines(unlist(qname), "TSS-TES.full-length.txt")

gene_matrices <- counts %>%
  group_by(gene_id) %>%
  group_split() %>%
  lapply(function(data) {
    data %>%
      dplyr::select(Promoter, PAS, read_counts) %>%
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
