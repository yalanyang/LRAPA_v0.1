library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("-b", "--bam_files"), type = "character", default = NULL, 
              help = "Comma-separated list of input BAM files", metavar = "character"),
  make_option(c("-p", "--pas"), type = "character", default = NULL, 
              help = "PAS reference TXT file with BED columns followed by annotations", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "output.txt", 
              help = "Output text file for the count matrix [default= %default]", metavar = "character")
)

## Rscript 3.PAS_count_Per_sample.R -b Neuron.filter.bam,Progenitor.filter.bam -p organoid.PAS.bed -o Neuron_Progentior.count.txt

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(dplyr)



# Check if required arguments are provided
if (is.null(opt$bam_files) || is.null(opt$pas)) {
  print_help(opt_parser)
  stop("Both --bam_files and --reference must be supplied", call. = FALSE)
}

# Function to read and process the PAS reference TXT file
read_pas_reference <- function(reference_file) {
  if (!file.exists(reference_file)) {
    stop("Reference file not found: ", reference_file)
  }
  # Read the TXT file into a data frame
  pas_df <- read.table(reference_file, header = TRUE, stringsAsFactors = FALSE)
  pas_df$Start=pas_df$Start+1
  # Create a GRanges object from the BED columns
  pas_gr <- GRanges(seqnames = pas_df$Chr,
                    ranges = IRanges(start = pas_df$Start, end = pas_df$End),
                    strand = pas_df$Strand)
  # Store additional annotation columns in metadata of GRanges object
  pas_gr$annotations <- pas_df[, c(10,11,9,14)]
  return(pas_gr)
}

# Read PAS reference TXT file
pas_reference <- read_pas_reference(opt$pas)

# Ensure pas_reference is a GRanges object
if (!inherits(pas_reference, "GRanges")) {
  stop("PAS reference file is not a valid GRanges object")
}

# Split the BAM files by comma
bam_files <- unlist(strsplit(opt$bam_files, ","))

# Function to extract PAS sites from BAM file
extract_pas_sites <- function(bam_file) {
  if (!file.exists(bam_file)) {
    stop("BAM file not found: ", bam_file)
  }
  bam <- readGAlignments(bam_file, param = ScanBamParam(what = c("rname", "strand", "pos", "qwidth")))
  
  gr <- GRanges(
    seqnames = seqnames(bam),
    ranges = IRanges(
      start = ifelse(strand(bam) == "+", end(bam), start(bam)),
      end = ifelse(strand(bam) == "+", end(bam), start(bam))
    ),
    strand = strand(bam)
  )
return(gr)
}

# Function to count PAS sites
count_pas_sites <- function(pas_gr, ref_gr) {
  overlaps <- findOverlaps(pas_gr,ref_gr)
  count_vector <- tabulate(subjectHits(overlaps), nbins = length(ref_gr))
  return(count_vector)
}

# Initialize count matrix dimensions
n_samples <- length(bam_files)
n_pas_sites <- length(pas_reference)
n_annotations <- length(colnames(pas_reference$annotations))

# Initialize count matrix with zeros
count_matrix <- matrix(0, nrow = n_pas_sites, ncol = n_samples)

# Process each BAM file
for (i in seq_along(bam_files)) {
  bam_file <- bam_files[i]
  pas_gr <- extract_pas_sites(bam_file)
  counts <- count_pas_sites(pas_gr,pas_reference)
  count_matrix[, i] <- counts
}

count_matrix <-cbind(count_matrix, pas_reference$annotations)
#rownames(count_matrix) <- paste0(seqnames(pas_reference), start(pas_reference), end(pas_reference), strand(pas_reference), sep = ":")
colnames(count_matrix) <- c(bam_files, colnames(pas_reference$annotations))

# Filter out rows with sum count = 0
sum_counts <- rowSums(count_matrix[, 1:n_samples])
count_matrix <- count_matrix[sum_counts > 0, ]

rownames(count_matrix) <- count_matrix$PAS_ID

# Write the count matrix to a text file
write.table(count_matrix, file = opt$output, sep = "\t", quote = FALSE, col.names = TRUE)

# Print a message indicating where the output was saved
cat("Count matrix saved to", opt$output, "\n")
