## update: 11-13-2024 
## To remove genes from a GTF file where the distance between two or more of its non-overlapping transcripts exceeds 10 kb
## and remove genes that have transcripts in different strands, such as Gm2004 in mouse refgene.

## update: 11-19-2024 
## To remove the first and last exon of each transcript

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
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
bam <- opt$input
pas.bed <- opt$pas
annot_path <-opt$gtf
output <- opt$output
#Examples: Rscript 5_PAS_exon_pair_V3.R -i Encode.TSS.TES2.bam -g /Users/yangyalan/results/Long-read-APA-pipeline/reference/hg38.refGene.gtf -p ../Encode.brain.new.PAS.bed -o TSS-exon.coordination.chiqtest.txt

library(GenomicRanges)
library(dplyr)
library(data.table)

library(stringr)

##prepare PAS sites
message("Peparing PAS sites")
pas <- fread(pas.bed, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pas <- pas[pas$Feature != "intergenic"]
pas <- pas %>% dplyr::filter(gene_count >= 10 & PAU >=0.05 & PAU <= 0.95) %>%  dplyr::select(Chr,Start,End,PAS_ID, gene_name,Strand) %>% arrange(Chr, Start, End)
pas$Start=pas$Start+1
pas.base <- tibble::as_tibble(pas) %>% dplyr::group_by(Chr,Strand,gene_name) %>% 
  dplyr::mutate(count = paste0(gene_name, ":T", sprintf("%02d", sequence(dplyr::n())))) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)


pas.anno <- DataFrame(count=elementMetadata(pas.base)$count,
                      PAS_ID_new=elementMetadata(pas.base)$PAS_ID)
colnames(pas.anno) <- c("tes_id","pas_id")
pas.anno <- as.data.frame(pas.anno)


##prepare alternative splcied exon reference
message("Peparing alternatively skipped exons")
refAnnotation <- rtracklayer::import.gff(annot_path)
valid_chromosomes <- paste0("chr", c(1:22, "X","Y"))
refAnnotation <- refAnnotation[seqnames(refAnnotation) %in% valid_chromosomes]
transcripts <- refAnnotation[refAnnotation$type == "transcript"]
gene_strands <- transcripts %>% 
  as.data.frame() %>%
  group_by(gene_name) %>%
  summarize(unique_strands = n_distinct(strand))

# Identify genes that appear on multiple strands
genes_on_multiple_strands <- gene_strands$gene_name[gene_strands$unique_strands > 1]

# Remove these genes from further analysis
transcripts <- transcripts[!transcripts$gene_name %in% genes_on_multiple_strands]

transcript_list <- split(transcripts, transcripts$gene_name)
# Identify genes where all distances between consecutive transcripts are <= 100 kb
genes_to_keep <- sapply(transcript_list, function(gene_transcripts) {
  gene_transcripts <- sort(gene_transcripts)
  distances <- start(gene_transcripts)[-1] - end(gene_transcripts)[-length(gene_transcripts)]
  all(distances <= 10000)
})

valid_genes <- names(genes_to_keep[genes_to_keep])
refAnnotation <- refAnnotation[refAnnotation$gene_name %in% valid_genes]
refExons <- refAnnotation[refAnnotation$type == "exon"]
utr <- refAnnotation[refAnnotation$type=="3UTR",]

##prepare bam file, extract read coordinate and exon coordinate
bamAlignments <- GenomicAlignments::readGAlignments(bam, use.names = TRUE)
read_junctions <- GenomicAlignments::junctions(bamAlignments, use.mcols = TRUE)
valid_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
alignments <- GenomicRanges::GRanges(bamAlignments)
alignments$name <- names(bamAlignments)
names(alignments) <- NULL


get_ref_exon <- function(exon_in_reference) {
    exons_df <- as.data.frame(exon_in_reference)
    exons_df <-  exons_df %>% arrange(gene_name, transcript_id, start)
    ##remove the first and last exon of each transcript. Only retain the alternative spliced exons
    exons_df_out1 <- exons_df %>% group_by(gene_name, transcript_id) %>% filter(row_number() == 1 | row_number() == n())
    
    exon_counts <- exons_df %>% group_by(gene_name, start, end) %>%  summarize(transcript_count = n_distinct(transcript_id))
    max_transcript_counts <- exons_df %>% group_by(gene_name) %>% summarize(max_transcripts = n_distinct(transcript_id))
    alternative_exons <- exon_counts %>%  left_join(max_transcript_counts, by = "gene_name") %>% filter(transcript_count < max_transcripts)
    alternative_exons_df <- alternative_exons %>% inner_join(exons_df, by = c("gene_name", "start", "end"))
    alternative_exons <- unique(alternative_exons_df %>% select(seqnames,start,end,gene_name,strand))
    
    exons_df_out1_bed <- makeGRangesFromDataFrame(exons_df_out1, 
                                            seqnames.field = "seqnames", 
                                            start.field = "start", 
                                            end.field = "end", 
                                            strand.field = "strand", 
                                            keep.extra.columns = TRUE)
    

   granges_bed <- makeGRangesFromDataFrame(alternative_exons, 
                                            seqnames.field = "seqnames", 
                                            start.field = "start", 
                                            end.field = "end", 
                                            strand.field = "strand", 
                                            keep.extra.columns = TRUE)
    overlaps <- findOverlaps(granges_bed, granges_bed)
    # Identify the indices of overlapping exons 
    overlap_indices <- queryHits(overlaps)[subjectHits(overlaps) != queryHits(overlaps)]
    # Remove overlapping exons
    non_overlapping_granges <- granges_bed[-unique(overlap_indices)] 

    overlaps2 <- findOverlaps(non_overlapping_granges, exons_df_out1_bed)
    overlap_indices2 <- queryHits(overlaps2)[subjectHits(overlaps2) != queryHits(overlaps2)]
    non_overlapping_exon2 <- non_overlapping_granges[-unique(overlap_indices2)] 
    non_overlapping_exon2 <- as.data.frame(non_overlapping_exon2)
    SE <- non_overlapping_exon2 %>% arrange(seqnames, start, end)
    SE.ref <- tibble::as_tibble(SE) %>% dplyr::group_by(gene_id) %>% 
    dplyr::mutate(count = paste0(gene_id, ":SE", sprintf("%02d", sequence(dplyr::n())))) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    overlaps <- findOverlaps(SE.ref, utr, maxgap = 10)
    SE.ref <- SE.ref[-queryHits(overlaps)]
    return(SE.ref)
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

readTESassignment <- function(endsAlignements, TESCoordinate.base) {
  tesDb <- TESCoordinate.base
  ovlps <- findOverlaps(endsAlignements , tesDb , maxgap = 24)
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

get_exon_count <- function(read_junctions,exon_in_reference){
  read_junctions <- read_junctions[seqnames(read_junctions) %in% valid_chromosomes]
  gr <- unlist(read_junctions)
  gr$read_id <- names(gr)
  names(gr) <- NULL
  df <- as.data.frame(gr)
  df_result <- df %>%
        group_by(read_id) %>%
              arrange(start) %>%
                    mutate(next_start = lead(start), next_end = lead(end)) %>%
                        filter(!is.na(next_start)) %>%
                             select(read_id, seqnames, start = end, end = next_start, strand)

  df_result_clean <- df_result %>% filter(start <= end)
  SEs <- GRanges(
            seqnames = df_result_clean$seqnames,
            ranges = IRanges(start = df_result_clean$start, end = df_result_clean$end),
            strand = df_result_clean$strand,
            read_id = df_result_clean$read_id
     )
  
  ovlps <- GenomicRanges::findOverlaps(SEs , exon_in_reference, maxgap = 50)
  SEIds <- SE.ref[subjectHits(ovlps),]$count
  SEStarts <- GenomicRanges::start(SE.ref[subjectHits(ovlps),])
  SEAlignemnts2 <- SEs[queryHits(ovlps),]
  SEAlignemnts2$SE_id <- SEIds
  SEAlignemnts2$SEStarts <- SEStarts
  SEAlignemnts2$dist2Assignment <-
  abs(GenomicRanges::start(SEAlignemnts2) - SEAlignemnts2$SEStarts)
  SEAlignemnts2 <- SEAlignemnts2[order(SEAlignemnts2$dist2Assignment, decreasing = FALSE),]
  SEAlignemnts2 <- SEAlignemnts2 %>% as.data.frame(.) %>% dplyr::select(read_id, SE_id)
  colnames(SEAlignemnts2) <- c("name","SE_id")
  pairsTested <-  left_join(endsAlignemnts, SEAlignemnts2, by = "name") 
  pairsTested <- pairsTested %>% dplyr::filter(gsub("\\:.*", "", tes_id) == gsub("\\:.*", "", SE_id)) 
  pairsTested <- pairsTested %>% mutate(gene_id =  gsub(":.*", "", pairsTested$tes_id))
  pairsTested <- pairsTested %>%
  mutate(pairs_id = paste0(
    gene_id,
    ":",
    gsub(".*:", "", pairsTested$tes_id),
    ":",
    gsub(".*:", "", pairsTested$SE_id)
  ))
  countedPairs <- pairsTested %>% group_by(pairs_id) %>% tally() %>%  dplyr::rename(in_in = n)
  countedPairsFinal <-
  left_join(
    countedPairs,
    pairsTested %>%
      dplyr::distinct(pairs_id, .keep_all = TRUE) %>%
      dplyr::select(gene_id, pairs_id, tes_id, SE_id),
    by = "pairs_id"
  )
  return(countedPairsFinal)
}

get_read_count <- function(alignmentsFile , exon_in_reference){
  ovlps2 <- GenomicRanges::findOverlaps(alignmentsFile , exon_in_reference)
  SEIds2 <- exon_in_reference[subjectHits(ovlps2),]$count
  SEStarts2 <- GenomicRanges::start(exon_in_reference[subjectHits(ovlps2),])
  SEAlignemnts3 <- alignmentsFile[queryHits(ovlps2),]
  SEAlignemnts3$SE_id <- SEIds2
  SEAlignemnts3$SEStarts <- SEStarts2
  SEAlignemnts3$dist2Assignment <-
  abs(GenomicRanges::start(SEAlignemnts3) - SEAlignemnts3$SEStarts)
  SEAlignemnts3<- SEAlignemnts3[order(SEAlignemnts3$dist2Assignment, decreasing = FALSE),]
  SEAlignemnts3 <- SEAlignemnts3 %>% as.data.frame(.) %>% dplyr::select(name, SE_id)
  colnames(SEAlignemnts3) <- c("name","SE_id")
  pairsTested2 <-  left_join(endsAlignemnts, SEAlignemnts3, by = "name") 
  pairsTested2 <- pairsTested2 %>% dplyr::filter(gsub("\\:.*", "", tes_id) == gsub("\\:.*", "", SE_id)) 

  pairsTested2 <- pairsTested2 %>% mutate(gene_id =  gsub(":.*", "", pairsTested2$tes_id))
  pairsTested2 <- pairsTested2 %>%
  mutate(pairs_id = paste0(
    gene_id,
    ":",
    gsub(".*:", "", pairsTested2$tes_id),
    ":",
    gsub(".*:", "", pairsTested2$SE_id)
  ))

  countedPairs <- pairsTested2 %>% group_by(pairs_id) %>% tally() %>%  dplyr::rename(all = n)
  return(countedPairs)
}


SE.ref <- get_ref_exon(refExons)

message("counting polyA reads")
endsAlignemnts <-  prepareForCountEnds(alignments, 1)
endsAlignemnts <- readTESassignment(endsAlignemnts, pas.base)
endsAlignemnts <- endsAlignemnts %>% as.data.frame(.) %>% dplyr::select(name, tes_id)

message("counting reads that contain a given exon")
counted_exon_Pairs <- get_exon_count(read_junctions, SE.ref)
message("counting reads that cover a given exon")
counted_read_Pairs <- get_read_count(alignments, SE.ref)

countedPairsFinal <-  left_join(counted_exon_Pairs, counted_read_Pairs, by = "pairs_id") 
countedPairsFinal <- countedPairsFinal %>% filter(all >= 10)
countedPairsFinal$in_out <- countedPairsFinal$all- countedPairsFinal$in_in
counts <- countedPairsFinal %>% group_by(gene_id) %>%  filter(n_distinct(tes_id) > 1)  %>%  ungroup()
counts <- counts %>% select(SE_id,tes_id,in_in,in_out)
counts <- counts %>% mutate(in_in = pmax(in_in, 0), in_out = pmax(in_out, 0))
counts <- counts %>%  group_by(SE_id) %>% filter(dplyr::n() > 1) %>% ungroup()
counts <-  left_join(counts, pas.anno, by = "tes_id") 


merged_data <- counts %>% select(SE_id,pas_id) %>%
  group_by(SE_id) %>%
  summarize(pas_id = str_c(pas_id, collapse = ";"))

write.table(counts,"PAS-exon.coordination.count.txt",sep="\t",quote=F, row.names = F)

chisq_test <- function(counts){
    gene_matrices <- counts %>%
        group_by(SE_id) %>%
        group_split() %>%
        lapply(function(data) {
        data %>%
        select(SE_id, tes_id, in_in, in_out,pas_id) %>%
         replace(is.na(.), 0)  # 将 NA 替换为 0
  })

names(gene_matrices) <- counts$SE_id %>% unique()

message("Peforming statistics analysis")
chisq_results <- lapply(gene_matrices, function(matrix) {
  matrix <- as.matrix(matrix[,c(3,4)])
  if (nrow(matrix) >= 2) {
    chisq.test(matrix)$p.value
  } else {
    NA  
  }
})

chisq_df <- data.frame(
         SE_id = names(chisq_results),
          p_value = unlist(chisq_results),
        
          stringsAsFactors = FALSE
   )
  chisq_df <- chisq_df[order(chisq_df$p_value), ]
  chisq_df$FDR <- p.adjust(chisq_df$p_value, method = "BH")
  chisq_df$SE_id <-  rownames(chisq_df)
  
  SE.ref2 <- as.data.frame(SE.ref)
  SE.ref2$exon <- paste0(SE.ref2$seqnames,":", SE.ref2$start,"-",SE.ref2$end, ":", SE.ref2$strand)
  SE.ref2$SE_id <-  SE.ref2$count
  SE.ref2 <- SE.ref2  %>%  select(SE_id,exon)
              
  results<- chisq_df %>% left_join(SE.ref2,by="SE_id")
  results<- results %>% left_join(merged_data,by="SE_id")
     return(results)
}

PAS_SE_co <- chisq_test(counts)
PAS_SE_co <- PAS_SE_co[complete.cases(PAS_SE_co[, c("p_value", "FDR")]), ]
PAS_SE_co$sig <- ifelse(PAS_SE_co$FDR < 0.05, "TRUE", "FALSE")
PAS_SE_co$gene_id <- sapply(strsplit(as.character(PAS_SE_co$SE_id), ":"), function(x) x[1])

write.table(PAS_SE_co, output, sep="\t",quote=F,row.names = F)

