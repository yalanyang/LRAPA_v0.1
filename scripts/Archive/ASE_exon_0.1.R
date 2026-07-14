library(optparse)
option_list <- list(
  make_option(c("-i", "--bam_file"), type = "character", default = NULL, 
              help = "input BAM files", metavar = "character"),
  make_option(c("-s", "--se"), type="character", default=NULL,
              help="Reference annotation gtf file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="ASE_exon.diff.txt",
              help="Output significant SEs", metavar="file"),
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
se.anno <-opt$se
output <- opt$output
bam_files <- unlist(strsplit(opt$bam_files, ","))

#Examples: Rscript ASE_exon_0.1.R -i B6/13.bam,SPRET/13.bam -s /Users/yangyalan/results/LRAPA_v0.1/references/SE/mouse.refSE_combine_gencode_refgene.txt -o S1_C13.ASE_exon.diff.txt

library(GenomicRanges)
library(dplyr)
library(data.table)
library(stringr)

message("Loading Skipped exons")
SE.anno<- fread(se.anno, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

SE.ref <- GRanges(
  seqnames = SE.anno$seqnames,
  ranges = IRanges(start = SE.anno$start, end = SE.anno$end),
  strand = SE.anno$strand,
  value.group_name = SE.anno$gene_name,
  count = SE.anno$count
)

SE.ref2 <- as.data.frame(SE.ref)
SE.ref2$exon <- paste0(SE.ref2$seqnames,":", SE.ref2$start,"-",SE.ref2$end, ":", SE.ref2$strand)
SE.ref2$SE_id <-  SE.ref2$count
SE.ref2 <- SE.ref2  %>%  select(SE_id,exon)


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
  countedSE <- SEAlignemnts2 %>% group_by(SE_id) %>% tally() %>%  dplyr::rename(inclusion = n)
  return(countedSE)
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
  
  countSE_gene <- SEAlignemnts3 %>% group_by(SE_id) %>% tally() %>%  dplyr::rename(total_counts = n)
  return(countSE_gene)
}



##prepare bam file, extract read coordinate and exon coordinate
countedFinal2 <- as.data.frame(0, nrow = 4 * length(SE.ref), ncol = length(bam_files)*2+1)

for (i in seq_along(bam_files)) {
  bam <- bam_files[i]
bamAlignments <- GenomicAlignments::readGAlignments(bam, use.names = TRUE)
read_junctions <- GenomicAlignments::junctions(bamAlignments, use.mcols = TRUE)
alignments <- GenomicRanges::GRanges(bamAlignments)
alignments$name <- names(bamAlignments)
names(alignments) <- NULL

message("counting reads that contain a given exon")
counted_SE <- get_exon_count(read_junctions, SE.ref)
message("counting reads that cover a given exon")
counted_SE_gene <- get_read_count(alignments, SE.ref)
countedFinal <-  left_join(counted_SE, counted_SE_gene, by = "SE_id")
#countedFinal <- countedPairsFinal %>% filter(all >= 10)
countedFinal$exclusion <- countedFinal$total_counts- countedFinal$inclusion
countedFinal$PSI <- countedFinal$inclusion / countedFinal$total_counts
colnames(countedFinal) <- c("SE_id", paste0(bam_files[i],"_inclusion"), paste0(bam_files[i],"_totalcounts"), paste0(bam_files[i],"_exclusion"), paste0(bam_files[i],"_PSI"))
SE.ref2 <- left_join(SE.ref2, countedFinal, by = "SE_id")
}

selected_cols <- grep("_all$", names(SE.ref2), value = TRUE)
countedFinal2 <- SE.ref2 %>%
  filter(if_all(all_of(selected_cols), ~ . > 10))

selected_cols <- grep("_exclusion$", names(countedFinal2), value = TRUE)
countedFinal2 <- countedFinal2 %>%
  filter(if_all(all_of(selected_cols), ~ . > 0))

countedFinal2$dPSI <- countedFinal2[,6]-countedFinal2[,10]
results <- list()

p_values <- lapply(1:nrow(countedFinal2), function(i) {
  matrix <- matrix(
    c(countedFinal2[i,3], countedFinal2[i,7],
      countedFinal2[i,5], countedFinal2[i,9]),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("B6", "Cast"), c("Inclusion", "Exclusion"))
  )
   chisq.test(matrix)$p.value
})

countedFinal2$pvalue <- p_values
countedFinal2$pvalue <- as.numeric(countedFinal2$pvalue)
countedFinal2 <- countedFinal2[order(countedFinal2$pvalue), ]
countedFinal2$adj_pvalue <- p.adjust(countedFinal2$pvalue, method = "BH")
countedFinal2$sig <- ifelse(abs(countedFinal2$dPSI) >= 0.1 & countedFinal2$adj_pvalue <= 0.05 , "TRUE", "FALSE")
countedFinal2$gene_id <- sapply(strsplit(as.character(countedFinal2$SE_id), ":"), function(x) x[1])

write.table(countedFinal2, "ASE_exon.diff.txt", sep = "\t", row.names = FALSE, quote = FALSE)

