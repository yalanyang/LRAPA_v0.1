## update: 11-13-2024 
## To remove genes from a GTF file where the distance between two or more of its non-overlapping transcripts exceeds 10 kb
## and remove genes that have transcripts in different strands, such as Gm2004 in mouse refgene.

##  update: 03-28-2025.V5
## For chisq.test,find the TSS-PAS pair with the highest contribution (dominant), returns the contribution and residual of dominant TSS-PAS pair.

##  update: 06-17-2026.V6
## using two PAs with most counts, annotation the type of coupling: ALE or TUTR

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
#Examples: Rscript 5_PAS_TSS_pair_V5.R -i Encode_adult/encode.mapping.unique.filter.bam -s ../prepare_TSS/Human.hg38.refTSS.combine.txt -p Encode_adult/Encode.brain.new.PAS.bed -o TSS-TES.encode.coordination.chiqtest.0617.txt
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
colnames(TSS.anno) <- c("tss_id","tss")

message("preparing poly(A) sites")
pas <- fread(pas.bed, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pas <- pas[pas$Feature != "intergenic"]
pas <- pas %>% dplyr::filter(gene_count >= 10 & PAU >=0.05 & PAU <= 0.95) %>%  dplyr::select(Chr,Start,End,PAS_ID, gene_name,Strand,UTR_id) %>% arrange(Chr, Start, End)
pas$Start=pas$Start+1
pas.base <- tibble::as_tibble(pas) %>% dplyr::group_by(Chr,Strand,gene_name) %>% 
  dplyr::mutate(count = paste0(gene_name, ":PA", sprintf("%02d", sequence(dplyr::n())))) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

pas.anno <- DataFrame(count=elementMetadata(pas.base)$count,
                      PAS_ID=elementMetadata(pas.base)$PAS_ID,
                      UTR_ID=elementMetadata(pas.base)$UTR_id)
write.table(pas.anno,"pas.anno.txt",sep="\t",quote=F, row.names = F)

colnames(pas.anno) <- c("pas_id","pas","utr_id")
pas.anno <- as.data.frame(pas.anno)
message("poly(A) database sucessfully created")

message("reading bam file")
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
  startsAlignemnts2$tss_id <- promoterIds
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
  endsAlignemnts2$pas_id <- tesIds
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
    dplyr::select(name, tss_id)
  # make single nt ends
  endsAlignemnts <-  prepareForCountEnds(alignmentsFile, 1)
  endsAlignemnts <- readTESassignment(endsAlignemnts, PAS)
  endsAlignemnts <- endsAlignemnts %>% as.data.frame(.) %>% dplyr::select(name, pas_id)
  # make pairs
  pairsTested <- left_join(endsAlignemnts, startAlignments, by = "name") 

  pairsTested <- pairsTested %>% dplyr::filter(gsub("\\:.*", "", pas_id) == gsub("\\:.*", "", tss_id)) 
  
  pairsTested <- pairsTested %>% mutate(gene_id =  gsub(":.*", "", pairsTested$tss_id))

  pairsTested <- pairsTested %>%
    mutate(pairs_id = paste0(
      gene_id,
      ":",
      gsub(".*:", "", pairsTested$tss_id),
      ":",
      gsub(".*:", "", pairsTested$pas_id)
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
counts <- counts %>% tidyr::separate(pairs_id, into = c("gene_id", "tss_id", "pas_id"), sep = ":") %>% filter(read_counts > 5)

gene_summary <- counts %>%
  group_by(gene_id) %>%
  summarize(
    num_tss = n_distinct(tss_id),
    num_PAS = n_distinct(pas_id)
  )

genes_with_multiple <- gene_summary %>%
  filter(num_tss > 1, num_PAS > 1) %>%
  pull(gene_id)

counts <- counts %>%
  filter(gene_id %in% genes_with_multiple) 


write.table(counts,"TSS-polyA.coordination.count.txt",sep="\t",quote=F, row.names = F)

counts_anno <- counts 
counts_anno$tss_id <- paste0(counts_anno$gene_id, ":", counts_anno$tss_id)
counts_anno$pas_id <- paste0(counts_anno$gene_id, ":", counts_anno$pas_id)
counts_anno <- left_join(counts_anno,TSS.anno, by="tss_id")
counts_anno <- left_join(counts_anno,pas.anno, by="pas_id")

merged_data1 <- counts_anno %>% dplyr::select(gene_id,tss) %>% unique() %>% group_by(gene_id) %>%
  summarize(tss = str_c(tss, collapse = ";"))

merged_data2 <- counts_anno %>% dplyr::select(gene_id,pas) %>% unique() %>% group_by(gene_id) %>%
  summarize(pas = str_c(pas, collapse = ";"))

qname <- countData$pairsTested$name
writeLines(unlist(qname), "TSS-TES.full-length.txt")

library(dplyr)
library(tidyr)

top2_pas <- counts_anno %>%
  group_by(gene_id, pas) %>%
  summarise(total_count = sum(read_counts, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(gene_id) %>%
  arrange(desc(total_count), .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  mutate(rank = paste0("PA", row_number())) %>%
  ungroup() %>%
  select(gene_id, rank, pas) %>%
  pivot_wider(
    names_from = rank,
    values_from = pas
  )

top2_pas <- left_join(top2_pas,pas.anno, by=c("PA1"="pas"))
top2_pas <- left_join(top2_pas,pas.anno, by=c("PA2"="pas"))
top2_pas <- top2_pas %>% dplyr::select(gene_id, PA1, PA2, utr_id.x, utr_id.y)

gene_matrices <- counts %>%
  group_by(gene_id) %>%
  group_split() %>%
  lapply(function(data) {
    data %>%
      dplyr::select(tss_id, pas_id, read_counts) %>%
      pivot_wider(names_from = tss_id, values_from = read_counts) %>%
      replace(is.na(.), 0)
  })
#names(gene_matrices) <- counts$gene_id %>% unique()
names(gene_matrices) <- counts %>%
  group_by(gene_id) %>%
  group_keys() %>%
  pull(gene_id)


get_top_contributor <- function(matrix) {
  PAS_names <- matrix[[1]]
  matrix <- as.matrix(matrix[, -1]) 
  rownames(matrix) <- PAS_names
  if (nrow(matrix) >= 2 && ncol(matrix) >= 2) {
    chisq <- chisq.test(matrix)
    # Compute contribution
    contrib <- (chisq$observed - chisq$expected)^2 / chisq$expected
    residuals <- chisq$residuals
    # Filter only positive residuals
    contrib[residuals <= 0] <- NA  # Remove non-positive residuals
    # Handle edge case: No positive residuals
    if (all(is.na(contrib))) {
      return(tibble(
        p_value = chisq$p.value,
        TSS = NA,
        PAS = NA,
        contribution = NA
      ))
    }
    
    # Get PAS (row) and TSS (column) names

    idx <- which(!is.na(contrib), arr.ind = TRUE)
    vals <- contrib[idx]
    ord <- order(vals, decreasing = TRUE)
    max_idx <- idx[ord[1], , drop = FALSE]

      result <- tibble(
      p_value = chisq$p.value,
      TSS = colnames(matrix)[max_idx[,"col"]],
      PAS = rownames(matrix)[max_idx[,"row"]],
      contribution = contrib[max_idx]
    )
    
    return(result)
  } else {
    return(tibble(
      p_value = NA,
      TSS = NA,
      PAS = NA,
      contribution = NA
    ))
  }
}

top_contributors <- lapply(gene_matrices, get_top_contributor)
top_contributors_df <- bind_rows(lapply(names(top_contributors), function(gene) {
  top_contributors[[gene]] %>% mutate(gene_id = gene)
}))


top_contributors_df <- top_contributors_df[, c("gene_id", "p_value", "TSS", "PAS", "contribution")]

top_contributors_df$tss_id <- paste0(top_contributors_df$gene_id,":",top_contributors_df$TSS)
top_contributors_df$pas_id <- paste0(top_contributors_df$gene_id,":",top_contributors_df$PAS)

top_contributors_df <- left_join(top_contributors_df,TSS.anno, by="tss_id")
top_contributors_df <- left_join(top_contributors_df,pas.anno, by="pas_id")
                 
top_contributors_df$dominant_pair_coor <- paste(top_contributors_df$tss, top_contributors_df$pas, sep = "|")
top_contributors_df$dominant_pair <- paste(top_contributors_df$TSS, top_contributors_df$PAS, sep = "-")

top_contributors_df <- left_join(top_contributors_df,top2_pas, by="gene_id")

top_contributors_df$annotation <- ifelse(
  top_contributors_df$utr_id.x == top_contributors_df$utr_id.y & !(top_contributors_df$utr_id.x == "." & top_contributors_df$utr_id.y == "."),
  "TUTR",
  "ALE"
)
top_contributors_df <- top_contributors_df[, c("gene_id", "p_value", "dominant_pair", "dominant_pair_coor", "contribution","PA1","PA2","annotation")]

chisq_df <- top_contributors_df[order(top_contributors_df$p_value), ]
chisq_df$FDR <- p.adjust(chisq_df$p_value, method = "BH")
chisq_df$sig <- ifelse(chisq_df$FDR < 0.05, "TRUE", "FALSE")
chisq_df <- left_join(chisq_df,merged_data1, by="gene_id")
chisq_df <- left_join(chisq_df,merged_data2, by="gene_id")

write.table(chisq_df, output, sep="\t",quote=F,row.names = F)




