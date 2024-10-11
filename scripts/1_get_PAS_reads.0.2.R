library(optparse)
#Rscript $codes/1_get_PAS_reads.0.2.R -i test.fl_flipped.unique.bam -r $reference/GRCh38.primary_assembly.genome.fa -o test.HQ.qname.txt

# Define command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input BAM file", metavar = "character"),
  make_option(c("-o", "--hq"), type = "character", default = "qname_HQ.txt", 
              help = "Output BAM file for high quality reads [default= %default]", metavar = "character"),
  make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="Reference genome FASTA file", metavar="file"),
  make_option(c("-f", "--flank_size"), type="integer", default=10, 
              help="Flank size on either side of the mapped end position for internal priming [default= %default]", metavar="integer"),
  make_option(c("-c", "--cores"), type = "integer", default = 4, 
              help = "Number of cores to use for parallel processing [default= %default]", metavar = "integer"),
  make_option(c("-b", "--batch_size"), type = "integer", default = 10000, 
              help = "Number of reads to process in each batch [default= %default]", metavar = "integer"),
  make_option(c("-a", "--adjacent_size"), type="integer", default=40, 
              help="Flank size to search for PAS singal [default= %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input BAM file must be supplied", call. = FALSE)
}



library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(optparse)
library(doParallel)
library(foreach)


adjacent_size <- opt$adjacent_size
reference_genome <- opt$reference
flank_size <- opt$flank_size

# Function to check if a sequence is mostly adenines (or thymines for reverse strand)
is_mostly_A <- function(sequence, reverse_strand = FALSE) {
  base <- ifelse(reverse_strand, "T", "A")
  sequence_length <- nchar(sequence)
  
  if (sequence_length <= 20) {
    base_count <- sum(strsplit(sequence, "")[[1]] == base)
    return(base_count / sequence_length >= 0.95)
  } 
  else {
    if (reverse_strand){
      first_20 <- substr(sequence, sequence_length-19, sequence_length)
      rest <- substr(sequence, sequence_length-min(sequence_length,40)+1, sequence_length-20)
    }
    else{
      first_20 <- substr(sequence, 1, 20)
      rest <- substr(sequence, 21, min(sequence_length,40))
    }
    
    first_20_count <- sum(strsplit(first_20, "")[[1]] == base)
    rest_count <- sum(strsplit(rest, "")[[1]] == base)
    
    return((first_20_count / 20 >= 0.8) && (rest_count / nchar(rest) >= 0.95))
  }
}


# Function to check for stretch of adenosines considering strand
stretch_size = 6
has_six_ad_adenosines <- function(sequence) {
  return(grepl(paste(rep("A", stretch_size), collapse = ""), sequence))
}



# Define PAS hexamers
pas_hexamers <- c("AATAAA", "TTTAAA", "AAGAAA", "AACAAA", "TATAAA", "AATGAA", 
                  "ATTAAA", "AGTAAA", "AATATA", "CATAAA", "ACTAAA", "GATAAA")

# Function to identify PAS hexamer in a sequence
identify_pas <- function(sequence) {
  for (hexamer in pas_hexamers) {
    if (grepl(hexamer, sequence)) {
      return(hexamer)
    }
  }
  return(NA)
}


# Load reference genome
reference <- readDNAStringSet(reference_genome)
# Adjust chromosome name if necessary
names(reference) <- sub(" +.*", "", names(reference))  # Extract the correct chromosome name


parse_cigar <- function(cigar_str) {
  lengths <- as.numeric(unlist(regmatches(cigar_str, gregexpr("\\d+", cigar_str))))
  operations <- unlist(regmatches(cigar_str, gregexpr("[A-Z]", cigar_str)))
  return(list(lengths = lengths, operations = operations))
}

# Function to process each read
process_read <- function(read) {
  cigartuples <- cigar(read)
  seq <- mcols(read)$seq
  qname <-  mcols(read)$qname
  strand <- as.character(strand(read)[1])
  reverse_strand <- (strand == "-")
  forward_strand <- (strand == "+")
  rname <- seqnames(read)[1]
  start <- start(read)
  end <- end(read)
  parsed_cigar <- parse_cigar(cigartuples)
  operations <- parsed_cigar$operations
  lengths <- parsed_cigar$lengths
    
  if (length(operations) > 0 && length(lengths) > 0 && as.vector(rname %in% names(reference)) && grepl("N", cigartuples)) {
    if (operations[1] == "S" && reverse_strand) {
      softclip_start <- lengths[1]
      seq1 <- substr(seq, 1, softclip_start)
      
      flanking_start <- max(start - flank_size, 1)
      flanking_end <- min(start + flank_size - 1, length(reference[[rname]]))
      seq2 <- as.character(reverseComplement(subseq(reference[[rname]], start=flanking_start, end=flanking_end)))
      
      flank_start <- start + 1
      flank_end <- min(start + adjacent_size,length(reference[[rname]]))
      seq3 <- as.character(reverseComplement(subseq(reference[[rname]], start=flank_start, end=flank_end)))
      if (is_mostly_A(seq1, reverse_strand) && !has_six_ad_adenosines(seq2) && !is.na(identify_pas(seq3))) {
        return(qname)
      }
    }
    
    else if (operations[length(operations)] == "S" &&  forward_strand) {
      softclip_end <- lengths[length(lengths)]
      seq1 <- substr(seq, nchar(seq) - softclip_end + 1, nchar(seq))
      
      flanking_start <- max(end - flank_size + 1, 1)
      flanking_end <- min(end + flank_size, length(reference[[rname]]))
      seq2 <- as.character(subseq(reference[[rname]], start=flanking_start, end=flanking_end))
      flank_start <- max(end - adjacent_size, 1)
      flank_end <- end - 1
      seq3 <- as.character(subseq(reference[[rname]], start=flank_start, end=flank_end))
      if (is_mostly_A(seq1, reverse_strand) && !has_six_ad_adenosines(seq2) && !is.na(identify_pas(seq3))) {
        return(qname)
    }
    }  
    #else {
    #  return("LQ")
    #}
  }
  #else{
   # return("LQ")
  #}
  }
  

registerDoParallel(cores = opt$cores)

# Process each read and categorize it in parallel
process_batch <- function(batch) {
 foreach(i = seq_along(batch), .combine = rbind, .packages = c("GenomicAlignments", "Biostrings")) %dopar% {
    qname <-process_read(batch[i])
    list(HQ = qname)
  }
}

# Open the BAM file for reading
bamIn <- BamFile(opt$input, yieldSize = opt$batch_size)
open(bamIn)

# Process reads in batches
HQ_reads <- list()
count=0
while (TRUE) {
 message(paste0("Processed read:", count))
  reads <- readGAlignments(bamIn, param = ScanBamParam(what = c("seq","qname")))
  if (length(reads) == 0) break
  batch_result <- process_batch(reads)
  HQ_reads <- list(HQ_reads, batch_result)
  count=count+10000
}

#HQ_bam <- do.call(c, HQ_reads)
#export(HQ_bam, BamFile(opt$hq, "wb"))
writeLines(unlist(HQ_reads), opt$hq)
# Print the number of reads without splicing junctions and hexamers
cat("Number of high quality reads :", nrow(unlist(HQ_reads)), "\n")


# Close parallel backend
stopImplicitCluster()

                        