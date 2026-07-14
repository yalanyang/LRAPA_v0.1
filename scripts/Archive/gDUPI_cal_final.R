library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("-c", "--count"), type = "character", default = NULL, 
              help = "Comma-separated list of input count file", metavar = "character"),
  make_option(c("-p", "--pas"), type = "character", default = NULL, 
              help = "PAS reference TXT file with BED columns followed by annotations", metavar = "character")
)

## Rscript gDUPI_cal_final.R -c fetal.brain.sc.count.txt -p fetal.brain.sc.PAS.bed 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if input files are provided
if (is.null(opt$count) || is.null(opt$pas)) {
  print_help(opt_parser)
  stop("Count file and polyA reference file must be provided.", call.=FALSE)
}

count <- opt$count
pas_reference <- opt$pas

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

count <- read.table(count, sep = "\t", header = T, stringsAsFactors = FALSE)
PAS<-read.table(pas_reference, header = T,sep="\t")
PAS <- PAS %>% dplyr::filter(PAU>=0.01 & PAU <= 0.99 & Count >=10)
count <- count %>% filter(PAS_ID %in% PAS$PAS_ID)
rownames(count) <- count$PAS_ID
methods <- c("gene","ALS","TUTR")
sample_id <- colnames(count[, !colnames(count) %in% c("PAS_ID", "UTR_id", "Feature", "gene_name")])
for (method in methods){
if (method=="gene") {
  counts=count[,c("gene_name", "PAS_ID", sample_id)]
  colnames(counts)[1] <- "gene_id"
  colnames(counts)[2] <- "feature_id"
  counts <- as.data.frame(counts)
}else if(method=="ALS"){
  count <- count[count$Feature=="3UTR",]
  counts=count[,c("gene_name", "UTR_id", sample_id)]
  colnames(counts)[1] <- "gene_id"
  colnames(counts)[2] <- "feature_id"
  counts <- counts %>%
  group_by(feature_id) %>%
  summarise(across(all_of(sample_id), sum), .groups = 'drop')
  counts$gene_id <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[1])
  counts <- as.data.frame(counts)
}else if(method=="TUTR"){
  count <- count[count$Feature=="3UTR",]
  counts=count[,c("UTR_id", "PAS_ID", sample_id)]
  colnames(counts)[1] <- "gene_id"
  colnames(counts)[2] <- "feature_id"
  counts <- as.data.frame(counts)
  
} else {
  message("not available")
}


calculate_gDPAU <- function(df) {
  n <- length(df)
  if (n == 1) return(NA)
  p <- df / sum(df)
  gDPAU <- sum((seq_len(n) - 1) * p) / (n - 1)
  return(gDPAU)
}


if(method=="gene" | method=="TUTR"){
counts$start <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[2])
counts$strand <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[3])
counts <- counts %>%
  arrange(
    gene_id,
    strand,
    ifelse(strand == "-", -as.integer(start), as.integer(start))
  )
unique_sites <- counts %>% select(gene_id, feature_id) %>% distinct()
gene_site_counts <- unique_sites %>% group_by(gene_id) %>%  summarize(Site_Count = n())
filtered_genes <- gene_site_counts %>% filter(Site_Count > 1) %>% pull(gene_id)
counts <- counts %>% filter(gene_id %in% filtered_genes)
count_data <- counts %>%  select(-feature_id, -start, -strand) %>%  pivot_longer(-gene_id, names_to = 'sample_id', values_to = 'Count')
count_data <- count_data %>% group_by(gene_id, sample_id) %>% mutate(PAU = Count / sum(Count)) %>% ungroup()

}else if (method=="ALS"){
  counts$order <- counts$feature_id
  counts <- counts %>%
    separate(order, into = c(NA, NA, NA, NA, "utr"), sep = ":", convert = TRUE)
  
  counts <- counts %>%
    group_by(gene_id) %>%
    mutate(utr = factor(utr, levels = paste0("u", sort(as.numeric(sub("u", "", utr)))))) %>%
    arrange(gene_id, utr) %>%
    ungroup()
  
  unique_sites <- counts %>% select(gene_id, feature_id) %>% distinct()
  gene_site_counts <- unique_sites %>% group_by(gene_id) %>%  summarize(Site_Count = n())
  filtered_genes <- gene_site_counts %>% filter(Site_Count > 1) %>% pull(gene_id)
  counts <- counts %>% filter(gene_id %in% filtered_genes)
  
  count_data <- counts %>%  select(-feature_id, -utr) %>%  pivot_longer(-gene_id, names_to = 'sample_id', values_to = 'Count')
  
  count_data <- count_data %>% group_by(gene_id, sample_id) %>% mutate(PAU = Count / sum(Count)) %>% ungroup()
  
}else {
  message("not available")
}

dpau_matrix <- count_data %>%
  group_by(gene_id, sample_id) %>%
  summarize(DPAU = calculate_gDPAU(PAU), .groups = 'drop') %>%
  pivot_wider(names_from = sample_id, values_from = DPAU)

dpau_matrix <- as.data.frame(dpau_matrix)
#dpau_matrix <- dpau_matrix[rowSums(dpau_matrix[, 1:8] < 0.95,na.rm = TRUE) > 0, ]

write.table(dpau_matrix, paste0(method, ".dpau_matrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

