library(optparse)
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "the count matrix file", metavar = "character"),
  make_option(c("-s", "--sample"), type = "character", default = NULL,
              help = "sample information", metavar = "character"),
  make_option(c("-c", "--group1"), type = "character", default = NULL,
              help = "Comma-separated list of samples in group 1", metavar = "character"),
  make_option(c("-n", "--group2"), type = "character", default = NULL,
              help = "Comma-separated list of samples in group 2", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
##Usage: Rscript 4.1_diff.DRIMSeq_0.2.R -i PRJDB15555.count.txt -s sample.txt -c C -n T

library(dplyr)
library(tidyr)
library(tibble)

library(DRIMSeq)
library(stringr)
count <- read.table(opt$input, sep = "\t", header = T, stringsAsFactors = FALSE)
sample <- read.table(opt$sample, sep = ",", header = T, stringsAsFactors = FALSE)

group1_sample <- opt$group1
group2_sample  <- opt$group2
samples <- sample %>% dplyr::filter(sample$Group==group1_sample | sample$Group==group2_sample)
sample_id <- samples$sample_id

methods <- c("all","between","within")
for (method in methods){
if (method=="all") {
  counts=count[,c("gene_name", "PAS_ID", sample_id)]
  colnames(counts)[1] <- "gene_id"
  colnames(counts)[2] <- "feature_id"
  merged_data <- counts %>% 
    select(gene_id,feature_id) %>%
    group_by(gene_id) %>%
    summarize(feature_id = str_c(feature_id, collapse = ";"))
  counts <- as.data.frame(counts)
}else if(method=="between"){
  count <- count[count$Feature=="3UTR",]
  counts=count[,c("gene_name", "UTR_id", sample_id)]
  colnames(counts)[1] <- "gene_id"
  colnames(counts)[2] <- "feature_id"
  counts <- counts %>%
  group_by(feature_id) %>%
  summarise(across(all_of(sample_id), sum), .groups = 'drop')
  counts$gene_id <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[1])
  merged_data <- counts %>% 
    select(gene_id,feature_id) %>%
    group_by(gene_id) %>%
    summarize(feature_id = str_c(feature_id, collapse = ";"))
  counts <- as.data.frame(counts)
}else if(method=="within"){
  count <- count[count$Feature=="3UTR",]
  counts=count[,c("UTR_id", "PAS_ID", sample_id)]
  colnames(counts)[1] <- "gene_id"
  colnames(counts)[2] <- "feature_id"
  counts <- as.data.frame(counts)
  merged_data <- counts %>% select(gene_id,feature_id) %>%
    group_by(gene_id) %>%
    summarize(feature_id = str_c(feature_id, collapse = ";"))
  
} else {
  message("not available")
}

#ref:https://bioconductor.org/packages/release/bioc/vignettes/DRIMSeq/inst/doc/DRIMSeq.pdf
d <- dmDSdata(counts=counts, samples=samples)
#plotData(d)
n <- nrow(samples)
n.small <- trunc(nrow(samples)/2)
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=2,
              min_samps_gene_expr=n, min_gene_expr=2)
design_full <- model.matrix(~ Group, data = samples(d))
design_full
set.seed(123)
d <- dmPrecision(d, design = design_full)
head(mean_expression(d), 3)
common_precision(d)
head(genewise_precision(d))
#plotPrecision(d)
#library(ggplot2)
#ggp <- plotPrecision(d)
#ggp + geom_point(size = 4)
d <- dmFit(d, design = design_full, verbose = 1)

testgroup <- paste0("Group",group2_sample)

d <- dmTest(d, coef=testgroup)
res <- DRIMSeq::results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

head(res)
##calcualte DPAU

calculate_gDPAU <- function(df) {
  n <- length(df)
  if (n == 1) return(NA)
  p <- df / sum(df)
  gDPAU <- sum((seq_len(n) - 1) * p) / (n - 1)
  return(gDPAU)
}

group1 <- sample[sample$Group==group1_sample,"sample_id"]
group2 <- sample[sample$Group==group2_sample,"sample_id"]


if(method=="all" | method=="within"){
counts$start <- sapply(strsplit(as.character(counts$feature_id), "_"), function(x) x[3])
counts$strand <- sapply(strsplit(as.character(counts$feature_id), "_"), function(x) x[4])
counts <- counts %>%
  arrange(gene_id, 
          case_when(
            strand == '-' ~ desc(as.integer(start)),  # Convert start to integer if needed
            TRUE ~ as.integer(start)  # Convert start to integer if needed
          )
  )
unique_sites <- counts %>% select(gene_id, feature_id) %>% distinct()
gene_site_counts <- unique_sites %>% group_by(gene_id) %>%  summarize(Site_Count = n())
filtered_genes <- gene_site_counts %>% filter(Site_Count > 1) %>% pull(gene_id)
counts <- counts %>% filter(gene_id %in% filtered_genes)
count_data <- counts %>%  select(-feature_id, -start, -strand) %>%  pivot_longer(-gene_id, names_to = 'sample_id', values_to = 'Count')
count_data <- count_data %>% group_by(gene_id, sample_id) %>% mutate(PAU = Count / sum(Count)) %>% ungroup()

}else if (method=="between"){
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

dpau_matrix$gDPAU1 <-  rowMeans(dpau_matrix[group1],na.rm = TRUE)
dpau_matrix$gDPAU2 <-  rowMeans(dpau_matrix[group2],na.rm = TRUE)

dpau_matrix$dgDPAU <- dpau_matrix$gDPAU1 -dpau_matrix$gDPAU2

if(method=="all" | method=="within"){
dpau_matrix$Note <- ifelse(dpau_matrix$dgDPAU > 0, "lengthen", "shorten")
}else if (method=="between"){
  dpau_matrix$Note <- ifelse(dpau_matrix$dgDPAU > 0, "distal", "proximal")

}
res <- res %>% left_join(dpau_matrix, by ="gene_id")
res <- res %>% left_join(merged_data, by ="gene_id")
res$sig <- ifelse(abs(res$dgDPAU) >= 0.1 & res$adj_pvalue <= 0.05 , "Yes", "No")
write.table(res, paste0(group1_sample, '_', group2_sample, '.', method, ".diff.DRIMSeq.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}
