## The code is used for merging overlapping 3'UTRs in the hg38 refGene gtf file.
## The overlapping 3'UTRs from the same genes are merged.

library(ggpubr)
library(data.table)
library(GenomicAlignments)
gr <- fread("reference/hg38.refGene.gtf", sep = "\t", header = FALSE, 
            col.names = c("Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"))

gr <- gr[!grepl("*random|*alt|*fix|chrUn*|chrM",gr$Chromosome),]
gr <- gr[grepl("3UTR",gr$Feature),]
# Function to get features from GTF file
features <- gr
features[, c("gene_id") := tstrsplit(Attributes, ';', fixed = TRUE)[1]]
features[, gene_id := sub('gene_id "', '', gene_id)]
features[, gene_id := sub('"', '', gene_id)]
#features[, gene_name := sub(' gene_name "', '', gene_name)]
#features[, gene_name := sub('"', '', gene_name)]
features <- features[, c('Chromosome', 'Start', 'End','Feature','gene_id', 'Strand')]
features <- unique(features)
gr <- GRanges(
  seqnames = features$Chromosome,
  ranges = IRanges(start = features$Start, end = features$End),
  strand = features$Strand,
  Feature = features$Feature,
  gene_id = features$gene_id
  #gene_name = features$gene_name
)

# Split GRanges object by UTR_start
gr_list <- split(gr, gr$gene_id)

# Merge overlapping regions within each gene
merge_by_gene <- function(gr) {
  reduced_gr <- reduce(gr, min.gapwidth=5)
  strand_ordered <- ifelse(strand(reduced_gr) == "+", start(reduced_gr), -start(reduced_gr))
  ordered_indices <- order(strand_ordered)
  reduced_gr <- reduced_gr[ordered_indices]
  ids <- paste0("u", seq_along(reduced_gr))
  mcols(reduced_gr) <- mcols(gr[1,])
  mcols(reduced_gr)$ID <- ids
  return(reduced_gr)
}

merged_gr_list <- lapply(gr_list, merge_by_gene)
merged_gr <- unlist(GRangesList(merged_gr_list))
merged_gr <- sortSeqlevels(merged_gr)
merged_gr <- sort(merged_gr)

# Create a data frame from the merged GRanges object
merged_data <- data.frame(
  Chromosome = seqnames(merged_gr),
  Start = start(merged_gr),
  End = end(merged_gr),
  Strand = strand(merged_gr),
  Feature = mcols(merged_gr)$Feature,
  gene_id = mcols(merged_gr)$gene_id,
  UTR_id = paste0(mcols(merged_gr)$gene_id, ":", seqnames(merged_gr), ":", start(merged_gr), ":", end(merged_gr), ":", mcols(merged_gr)$ID)
)

write.table(merged_data, "hg38.refGene.3utr_merge.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=T)

merged_data$Start <- merged_data$Start-1
merged_data$score <- "."
merged_data <- merged_data[,c(1:3,5,7,4,6)]

write.table(merged_data, "reference/hg38.refGene.3utr_merge.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=F)


