## update: 11-13-2024 
## To remove genes from a GTF file where the distance between two or more of its non-overlapping transcripts exceeds 100 kb
## and remove genes that have transcripts in different strands, such as Gm2004 in mouse refgene.

## update: 11-19-2024 
## To remove the first and last exon of each transcript

library(optparse)
library(GenomicRanges)
library(dplyr)
library(data.table)
library(stringr)

option_list <- list(
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="Reference annotation gtf file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="Encode.V40.SE.ref.txt",
              help="Output SE reference based on gtf file", metavar="file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
annot_path <-opt$gtf
output <- opt$output

#Examples: Rscript prepare_SE.v1.R -g gencode.v40.annotation_utr.sorted.gtf -o Human.encode.V40.SE.ref.txt
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
  all(distances <= 100000)
})

valid_genes <- names(genes_to_keep[genes_to_keep])
refAnnotation <- refAnnotation[refAnnotation$gene_name %in% valid_genes]
refExons <- refAnnotation[refAnnotation$type == "exon"]
utr <- refAnnotation[refAnnotation$type %in% c("3UTR", "three_prime_utr"), ]

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
    #overlaps <- findOverlaps(granges_bed, granges_bed)
    # Identify the indices of overlapping exons 
    #overlap_indices <- queryHits(overlaps)[subjectHits(overlaps) != queryHits(overlaps)]
    # Remove overlapping exons
    #non_overlapping_granges <- granges_bed[-unique(overlap_indices)] 
    overlaps2 <- findOverlaps(granges_bed, exons_df_out1_bed)
    overlap_indices2 <- queryHits(overlaps2)[subjectHits(overlaps2) != queryHits(overlaps2)]
    non_overlapping_exon2 <- granges_bed[-unique(overlap_indices2)] 
    non_overlapping_exon2 <- as.data.frame(non_overlapping_exon2)
    SE <- non_overlapping_exon2 %>% arrange(seqnames, start, end)
    SE.ref <- tibble::as_tibble(SE) %>% dplyr::group_by(gene_name) %>% 
    dplyr::mutate(count = paste0(gene_name, ":SE", sprintf("%02d", sequence(dplyr::n())))) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    overlaps <- findOverlaps(SE.ref, utr, maxgap = 10)
    SE.ref <- SE.ref[-queryHits(overlaps)]
    return(SE.ref)
}


SE.ref <- get_ref_exon(refExons)
SE.ref <- data.frame(SE.ref)
write.table(SE.ref, output, sep="\t",quote=F,row.names = F)

