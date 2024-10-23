# Command line options
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input BAM file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output PAS file", metavar="file"),
  make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="Reference genome FASTA file", metavar="file"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="Reference annotation gtf file", metavar="file"),
  make_option(c("-u", "--utr"), type="character", default=NULL,
              help="Reference annotation 3utr file", metavar="file"),
  make_option(c("-f", "--flank_size"), type="integer", default=40, 
              help="Flank size to search for PAS [default= %default]", metavar="integer")
)

#Rscript 2_PAS_identification_Annotation_0.3.R -i organoid.HQ.sort.bam -r reference/GRCh38.primary_assembly.genome.fa -g reference/hg38.refGene.gtf -u reference/hg38.3utr_merge.bed -o organoid.PAS.bed


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if input files are provided
if (is.null(opt$input) || is.null(opt$reference)) {
  print_help(opt_parser)
  stop("Input BAM file and reference genome FASTA file must be provided.", call.=FALSE)
}

library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(Rsamtools)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(bedtoolsr)
library(plyranges)



# Assign command line arguments to variables
input_bam <- opt$input
reference_genome <- opt$reference
output <- opt$output
gtf_file <- opt$gtf
utr_file <- opt$utr
flank_size <- opt$flank_size

# Read BAM file into GAlignments object
reads <- readGAlignments(input_bam, param = ScanBamParam(what = c("seq", "strand", "cigar", "qname", "rname", "pos")))
utr_ref <- read.table(utr_file, header = F, stringsAsFactors = FALSE)

# Extract read end positions considering strand
polyA_sites <- GRanges(
  seqnames = seqnames(reads),
  ranges = IRanges(
    start = ifelse(strand(reads) == "+", end(reads), start(reads)),
    end = ifelse(strand(reads) == "+", end(reads), start(reads))
  ),
  strand = strand(reads)
)

count <- as.data.frame(polyA_sites) %>%
  group_by(seqnames, start, end, strand) %>%
  summarise(count = n(),.groups="drop") 

polyA_sites <- makeGRangesFromDataFrame(count, keep.extra.columns = TRUE)


#reference_genome <- "reference/GRCh38.primary_assembly.genome.fa"
reference <- readDNAStringSet(reference_genome)
names(reference) <- sub(" +.*", "", names(reference))  # Extract the correct chromosome name

sequences <- c()
for (i in seq_along(polyA_sites)) {
  chr <- as.character(seqnames(polyA_sites)[i])
  start <- start(polyA_sites)[i]
  end <- end(polyA_sites)[i]
  strand <- as.character(strand(polyA_sites)[i])
  if (chr %in% names(reference)) {
    seq <- subseq(reference[[chr]], start=start, end=end)
    if (strand == "-") {
      seq <- reverseComplement(seq)
    }
    sequences[i] <- as.character(seq)
  } else {
    # If chromosome is not found, return NA
    sequences[i] <- NA
  }
}
polyA_sites$seq <- sequences

polyA_sites <- sort(polyA_sites)
merged_gr <- reduce(polyA_sites, min.gapwidth = 24)

total_counts <- sapply(seq_along(merged_gr), function(i) {
  overlaps <- findOverlaps(polyA_sites, merged_gr[i])
  sum(polyA_sites$count[queryHits(overlaps)])
})


##select cleavage sites. 
selected_coords <- sapply(seq_along(merged_gr), function(i) {
  overlaps <- findOverlaps(polyA_sites, merged_gr[i])
  overlapping_ranges <- polyA_sites[queryHits(overlaps)]
  
  a_ranges <- overlapping_ranges[overlapping_ranges$seq == "A"]
  if (length(a_ranges) > 0 &&  max(a_ranges$count) >=10) {
    selected <- a_ranges[which.max(a_ranges$count)]
  } else {
    selected <- overlapping_ranges[which.max(overlapping_ranges$count)]
  }
  return(selected)
})

selected_coords <- do.call(c, selected_coords)

merged_gr$PAS_id <- paste0(seqnames(selected_coords), ":", 
                           end(selected_coords), ":", 
                           strand(selected_coords))
merged_gr$position <- end(selected_coords)
merged_gr$total_count <- total_counts

##only contained polyA_sites that are supported by at least two sites
clustered_sites <- merged_gr %>% filter(total_count>1)
message(paste0("PAS cluster done, total number of PAS clusters: ", length(clustered_sites$total_count)))


##hexamers annotation
# Define PAS hexamers
pas_hexamers <- c("AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AAGAAA","AATATA", 
                  "TTTAAA","CATAAA","AACAAA", "AATGAA", "ACTAAA", "GATAAA")

# Function to identify PAS hexamer in a sequence
identify_pas <- function(sequence) {
  for (hexamer in pas_hexamers) {
    if (grepl(hexamer, sequence)) {
      return(hexamer)
    }
  }
  return(NA)
}

pas_hexamers_for_sites <- list()

for (i in seq_along(clustered_sites)) {
  chr <- as.character(seqnames(clustered_sites[i]))
  end_pos <- clustered_sites$position[i]
  if (chr %in% names(reference)) {
    if (as.character(strand(clustered_sites[i]))=="+")
    {
      flank_start <- max(end_pos - flank_size, 1)
      flank_end <- end_pos - 1
      sequence <- as.character(subseq(reference[[chr]], start=flank_start, end=flank_end))
    }
    else
    {
      flank_start <- end_pos + 1
      flank_end <- min(end_pos + flank_size,length(reference[[chr]]))
      sequence <- as.character(reverseComplement(subseq(reference[[chr]], start=flank_start, end=flank_end)))
    }
  }
  
  pas_hexamer <- identify_pas(sequence)
#  if (!is.na(pas_hexamer)) {
    pas_hexamers_for_sites <- c(pas_hexamers_for_sites, pas_hexamer)
#  }
}

mcols(clustered_sites)$hexamer <- pas_hexamers_for_sites
clustered_sites <- clustered_sites %>% filter(!is.na(hexamer))

# Write clusters to a BED file
bed_data <- data.frame(
  chr = as.character(seqnames(clustered_sites)),
  start = start(clustered_sites) - 1,  # BED format uses 0-based start
  end = end(clustered_sites),
  hexamer=as.character(clustered_sites$hexamer),  # Optional name field for BED format
  PAS_id = clustered_sites$PAS_id,   # Optional score field for BED format
  strand = as.character(strand(clustered_sites)),
  count =clustered_sites$total_count
)

# Helper function to read and process GTF file
gr <- fread(gtf_file, sep = "\t", header = FALSE, 
              col.names = c("Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"))
gr <- gr[!grepl("*CDS",gr$Feature),]

feature_order <- c('3UTR',"stop_codon", "5UTR","start_codon",'exon','transcript',"intergenic")


# Function to get features from GTF file

  features <- gr[gr$Feature %in% feature_order, ]
  features[, "gene_id" := tstrsplit(Attributes, ';', fixed = TRUE)[1]]
  features[, gene_id := sub('gene_id "', '', gene_id)]
  features[, gene_id := sub('"', '', gene_id)]
  features$Start <- features$Start -1
  features <- features[, c('Chromosome', 'Start', 'End','Feature','gene_id','Strand')]


# Function to annotate BED file with GTF features
annotate_bed <- function(bed_file) {
  annotated_bed <- bt.intersect(a=bed_file, b=utr_ref, wa=T,wb=T,s=T,loj =T)

  annotated_bed <- annotated_bed[, c(1:7,12,14)]
  annotated_bed <- annotated_bed %>%
    arrange(V1, V2, V3) %>%
    distinct(V1, V2, V3, V6, .keep_all = TRUE)
  annotated_bed <- bt.intersect(a=annotated_bed, b=features,wa=T,wb=T,s=T,loj=T)
  annotated_bed <- annotated_bed[, c(1:8,13,14)]
  annotated_bed <- annotated_bed[!duplicated(annotated_bed), ]
  # Replace NA values and aggregate annotations by region

  annotated_bed$V13 <- ifelse(annotated_bed$V13==".", 'intergenic', annotated_bed$V13)
  
  
  colnames(annotated_bed) <- c('Chr', 'Start', 'End', "Hexamer", "PAS_id", 'Strand', "Count","UTR_id", 'Feature','gene_name')
  feature_priority <- setNames(seq_along(feature_order), feature_order)
  
  annotated_bed$feature_priority <- feature_priority[match(annotated_bed$Feature, feature_order)]
  
  annotated_bed$Feature[is.na(annotated_bed$Feature)] <- max(feature_priority) + 1
  
  unique_annotated_bed <- annotated_bed %>%
    arrange(Chr, Start, feature_priority) %>%
    distinct(Chr, Start, End, .keep_all = TRUE) %>% dplyr::select(-feature_priority)
  
  unique_annotated_bed$Feature <- ifelse(unique_annotated_bed$Feature=="transcript", 'intron', unique_annotated_bed$Feature)
  unique_annotated_bed$gene_name <- ifelse(unique_annotated_bed$Feature=="intergenic", paste0("intergenic",":", unique_annotated_bed$Chr,":", unique_annotated_bed$End), unique_annotated_bed$gene_name)

  
  Count_total <- aggregate(unique_annotated_bed$Count,by=list(type=unique_annotated_bed$gene_name),sum)
  colnames(Count_total) <- c("gene_name","gene_count")
  unique_annotated_bed = unique_annotated_bed %>% left_join(Count_total, by="gene_name")
  unique_annotated_bed$PAU <- unique_annotated_bed$Count/unique_annotated_bed$gene_count
  
  
  refGene <- gr[!grepl("*random|*alt|*fix|chrUn*|chrM",gr$Chromosome),]
  refGene <- refGene[grepl("transcript",refGene$Feature),]
  refGene$PAS <- ifelse(refGene$Strand == "+", refGene$End, refGene$Start)
  refGene$Start <- refGene$PAS-1
  refGene <- refGene %>% select(Chromosome,Start,PAS,Score,Frame,Strand)
  colnames(refGene) <- c('chrom', 'start', 'end', 'name', 'score', 'strand')
  refGene <- unique(refGene)
  refGene$name <- paste0(refGene$chrom, ":", refGene$end, ":", refGene$strand)
  refGene$start <- refGene$start-20
  refGene$end <- refGene$end+20
  unique_annotated_bed <- unique(bt.intersect(a=unique_annotated_bed, b=refGene, wa=T, wb=T, loj=T, s=T))
  unique_annotated_bed <- unique(unique_annotated_bed[,c(1:12,16)])
  #unique_annotated_bed$PAS <- paste0(unique_annotated_bed$V1, ":", unique_annotated_bed$V3, ":", unique_annotated_bed$V6)
  unique_annotated_bed$anno <- ifelse(unique_annotated_bed$V16 == ".", "not_RefGene", "yes_RefGene")
  unique_annotated_bed$PAS_ID_new <- ifelse(unique_annotated_bed$V16 == ".", unique_annotated_bed$V5, unique_annotated_bed$V16)
  unique_annotated_bed <- unique_annotated_bed %>% group_by(V1, V2, V3) %>%  slice(1) %>% ungroup()
  colnames(unique_annotated_bed) = c("Chr",	"Start",	"End",	"Hexamer",	"PAS_ID",	"Strand",	"Count","UTR_id",	"Feature",
                            "gene_name",	"gene_count",	"PAU", "annotated_PAS", "anno","annotated_PAS")
  
  return(unique_annotated_bed)
}

#Annotate regions
annotated_bed <- annotate_bed(bed_data)
# write the annotated BED file
#write.table(annotated_bed,'raw.pas.txt', quote=FALSE, sep="\t", row.names=FALSE, col.names=T)
write.table(annotated_bed, output, quote=FALSE, sep="\t", row.names=FALSE, col.names=T)

