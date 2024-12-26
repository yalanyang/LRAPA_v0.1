## 12-24-2024 Prepare TSS reference for coupling_analysis


library(optparse)
library(AnnotationDbi)
library(dplyr)
library(data.table)
library(tidyr)
library(GenomicRanges)
library(stringr)

option_list <- list(
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="Reference annotation gtf file", metavar="file"),
  make_option(c("-b", "--bin"), type = "integer", default = 50, 
              help = "tss window [default= %default]", metavar = "integer"),
  make_option(c("-o", "--output"), type="character", default="Encode.V40.TSS.ref.txt",
              help="Output TSS reference based on gtf file", metavar="file")
)

#Examples: Rscript prepare_TSS.v1.R -g /Users/yangyalan/results/Long-read-APA-pipeline/reference/gencode.v40.annotation_utr.sorted.gtf -o /Users/yangyalan/results/Long-read-APA-pipeline/reference/Human.encode.V40.TSS.ref.txt

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
annot_path <- opt$gtf
output <- opt$output

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
  
  annotation$gene_id_clean <- gsub("\\..*", "", annotation$gene_id)
  e2$gene_id_clean <- gsub("\\..*", "", e2$gene_id)
  
  gene_mapping <- data.frame(
    gene_id_clean = annotation$gene_id_clean,
    gene_name = annotation$gene_name
  )
  
gene_mapping <- unique(gene_mapping)
e2$gene_name <- gene_mapping$gene_name[match(e2$gene_id_clean, gene_mapping$gene_id_clean)]

  names(e2) <- NULL
  mcols(e2) <- mcols(e2)[, c(
    "exon_id", "exon_rank",
    "transcript_id", "gene_name", "type"
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
          GenomicRanges::split(tss.bins, ~gene_name),
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
  all(distances <= 100000)
})

valid_genes <- names(genes_to_keep[genes_to_keep])
refAnnotation <- refAnnotation[refAnnotation$gene_id %in% valid_genes]

refExons <- refAnnotation[refAnnotation$type == "exon"]
TSS_database <- prepareTSSDatabase(refExons, tss.window=50)
TSS.anno <- as.data.frame(TSS_database)
write.table(TSS.anno,output,sep="\t",quote=F, row.names = F)


