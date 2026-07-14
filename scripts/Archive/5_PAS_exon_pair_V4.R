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
  make_option(c("-s", "--se"), type="character", default=NULL,
              help="Reference annotation gtf file", metavar="file"),
  make_option(c("-p", "--pas"), type="character", default=NULL,
              help="Reference polyA.bed", metavar="file")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
bam <- opt$input
pas.bed <- opt$pas
se.anno <-opt$se
output <- opt$output

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

message("Loading Skipped exons")
SE.anno<- fread(se.anno, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
SE.anno <- SE.anno[SE.anno$gene_name %in% unique(pas$gene_name), ]

SE.ref <- GRanges(
  seqnames = SE.anno$seqnames,
  ranges = IRanges(start = SE.anno$start, end = SE.anno$end),
  strand = SE.anno$strand,
  value.group_name = SE.anno$gene_name,
  count = SE.anno$count
)

SE.anno$exon <- paste0(SE.anno$seqnames,":", SE.anno$start,"-",SE.anno$end, ":", SE.anno$strand)
SE.anno$SE_id <-  SE.anno$count
SE.anno <- SE.anno  %>%  select(SE_id,exon)

##prepare bam file, extract read coordinate and exon coordinate
message("Loading BAM file")
bamAlignments <- GenomicAlignments::readGAlignments(bam, use.names = TRUE)
read_junctions <- GenomicAlignments::junctions(bamAlignments, use.mcols = TRUE)
valid_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
bamAlignments <- bamAlignments[seqnames(bamAlignments) %in% valid_chromosomes]
read_junctions <- read_junctions[seqnames(read_junctions) %in% valid_chromosomes]
alignments <- GenomicRanges::GRanges(bamAlignments)
alignments$name <- names(bamAlignments)
names(alignments) <- NULL

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

get_exon_count <- function(read_junctions,exon_in_reference){
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
        
  results<- chisq_df %>% left_join(SE.anno,by="SE_id")
  results<- results %>% left_join(merged_data,by="SE_id")
     return(results)
}

PAS_SE_co <- chisq_test(counts)
PAS_SE_co <- PAS_SE_co[complete.cases(PAS_SE_co[, c("p_value", "FDR")]), ]
PAS_SE_co$sig <- ifelse(PAS_SE_co$FDR < 0.05, "TRUE", "FALSE")
PAS_SE_co$gene_id <- sapply(strsplit(as.character(PAS_SE_co$SE_id), ":"), function(x) x[1])

write.table(PAS_SE_co, output, sep="\t",quote=F,row.names = F)

