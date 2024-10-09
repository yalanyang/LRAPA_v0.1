library(optparse)
# Define command line arguments
option_list <- list(
  make_option(c("-b", "--bam_file"), type = "character", default = NULL, 
              help = "input BAM file", metavar = "character"),
  make_option(c("-p", "--pas"), type = "character", default = NULL, 
              help = "PAS reference TXT file with BED columns followed by annotations", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "count", 
              help = "Output folder for the count matrix [default= %default]", metavar = "character"),
  make_option(c("-w", "--barcode"), type = "character", default = NULL, 
              help = "whitelist barcodes", metavar = "character")
)

##Rscript 3_PAS_count_sc_sample_0.3.R -b RUN3.HQ.bam -p organiods.PAS.bed -o RUN3 -w barcode.RUN3.rev.txt

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$bam_file) || is.null(opt$pas) || is.null(opt$barcode)) {
  print_help(opt_parser)
  stop("Both --bam_file and --reference and --barcode must be supplied", call. = FALSE)
}

library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(dplyr)
library(tidyr)


barcode <- read.table(opt$barcode, header = T, stringsAsFactors = FALSE)
barcode_list <- barcode$barcode

# Function to read and process the PAS reference TXT file
read_pas_reference <- function(reference_file) {
  pas_df <- read.table(reference_file, header = TRUE, stringsAsFactors = FALSE)
  pas_df <- pas_df %>% dplyr::filter(Chr %in% paste0("chr", c(1:22, "X", "Y")))
  #pas_df <- pas_df %>% dplyr::filter(gene_count >= 10)
  pas_df <- pas_df %>% dplyr::filter(Count >= 5)
  pas_df <- pas_df %>% dplyr::filter(PAU > 0.01)
  
  pas_df$Start=pas_df$Start+1
  # Create a GRanges object from the BED columns
  pas_ref <- GRanges(seqnames = pas_df$Chr,
                    ranges = IRanges(start = pas_df$Start, end = pas_df$End),
                    strand = pas_df$Strand)
  # Store additional annotation columns in metadata of GRanges object
  pas_ref$annotations <- pas_df %>%  dplyr::select(gene_name,PAS_ID_new,Feature)
  return(pas_ref)
}

pas_reference <- read_pas_reference(opt$pas)


# Split the BAM files by comma
bam <- opt$bam_file

# Function to extract PAS sites from BAM file
extract_pas_sites <- function(bam_file) {
  bam <- readGAlignments(bam_file, param = ScanBamParam(what = c("qname"), tag = 'CB'))
  bam <- bam[mcols(bam)$CB %in% barcode_list]
  gr <- GRanges(
    seqnames = seqnames(bam),
    ranges = IRanges(
      start = ifelse(strand(bam) == "+", end(bam), start(bam)),
      end = ifelse(strand(bam) == "+", end(bam), start(bam))
    ),
    strand = strand(bam),
    qname = mcols(bam)$qname,
    CB = mcols(bam)$CB
  )
return(gr)
}

# Function to count PAS sites
count_pas_sites <- function(pas_reads, pas_ref) {
  overlaps <- findOverlaps(pas_reads,pas_ref)
  # Get all unique barcodes
  unique_barcodes <- unique(mcols(pas_reads)$CB)
  
  # Initialize a matrix with rows as PAS sites and columns as barcodes
  count_matrix <- matrix(0, nrow = length(pas_ref), ncol = length(unique_barcodes),
                         dimnames = list(seq_along(pas_ref), unique_barcodes))
  
  # Loop over each reference PAS site
  for (i in seq_along(pas_ref)) {
    # Get the PAS indices overlapping the current reference GRanges
    pas_idx <- queryHits(overlaps)[subjectHits(overlaps) == i]
    
    # If there are any overlaps, count the CBs and fill the matrix
    if (length(pas_idx) > 0) {
      cb_counts <- table(mcols(pas_reads)$CB[pas_idx])
      count_matrix[i, names(cb_counts)] <- as.integer(cb_counts)
    }
  }
  ID <- paste0(mcols(pas_ref)$annotations$gene_name,"_",mcols(pas_ref)$annotations$PAS_ID_new,"_", mcols(pas_ref)$annotations$Feature)
  rownames(count_matrix) <- ID

  return(count_matrix)
}

pas_reads <- extract_pas_sites(bam)
counts <- count_pas_sites(pas_reads, pas_reference)
cnt.out.path <- opt$output 
if (!dir.exists(cnt.out.path)) dir.create(cnt.out.path)
write.table(counts, file = paste0(cnt.out.path, "/count.matrix.txt"), sep = "\t", quote = FALSE, col.names = TRUE)

#sum_counts <- rowSums(counts[, 1:ncol(counts)])
#counts <- counts[sum_counts > 0, ]

#sum_counts_col <- colSums(counts[1:nrow(counts),])

#count_matrix <- counts[,sum_counts_col > 0]
count_matrix <- Matrix::Matrix(counts,sparse = TRUE)

#rownames(counts) <- counts$PAS_ID
pas <-  rownames(count_matrix)
gene <- mcols(pas_reference)$annotations$gene_name

feature <- data.frame(gene = gene, pas = pas)






Matrix::writeMM(count_matrix, file = paste0(cnt.out.path, "/matrix.mtx"))
utils::write.table(colnames(count_matrix), file = paste0(cnt.out.path, "/barcodes.tsv"),  quote = FALSE, row.names = FALSE, col.names = FALSE)
utils::write.table(feature, file = paste0(cnt.out.path, "/features.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)


# Print a message indicating where the output was saved
cat("Count matrix saved to", opt$output, "\n")


