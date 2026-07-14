library(optparse)
library(dplyr)
library(tidyr)
library(tibble)
library(DRIMSeq)
library(stringr)
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
##Usage: Rscript 4.1_diff.DRIMSeq_0.3.R -i PAS.count.txt -s sample.txt -c case -n control


count <- read.table(opt$input, sep = "\t", header = T, stringsAsFactors = FALSE)
sample <- read.table(opt$sample, sep = ",", header = T, stringsAsFactors = FALSE)

group1_sample <- opt$group1
group2_sample  <- opt$group2
samples <- sample %>% dplyr::filter(sample$Group==group1_sample | sample$Group==group2_sample)
sample_id <- samples$sample_id

##calcualte DPAU
calculate_gDPAU <- function(df) {
  n <- length(df)
  if (n == 1) return(NA)
  p <- df / sum(df)
  gDPAU <- sum((seq_len(n) - 1) * p) / (n - 1)
  return(gDPAU)
}

##calcualte WARM
calculate_WARM <- function(df) {
  WARM <- sum(df)
  return(WARM)
}


##calcualte MPRO
calculate_MPRO <- function(counts) {
  counts_MPRO <- counts
  counts_MPRO$group1 <- rowSums(counts_MPRO[, group1], na.rm = TRUE)
  counts_MPRO$group2 <- rowSums(counts_MPRO[, group2], na.rm = TRUE)
  counts_MPRO <- counts_MPRO %>%
    group_by(gene_id) %>%
    mutate(PAU1 = group1 / sum(group1, na.rm = TRUE), 
           PAU2 = group2 / sum(group2, na.rm = TRUE)) %>%
     ungroup()
  
  ## dPAU
  counts_MPRO$dPAU <- counts_MPRO$PAU1 - counts_MPRO$PAU2
  ## MPRO
  MPRO <- counts_MPRO %>%
    group_by(gene_id) %>%
    reframe({
      x <- dPAU
      n <- length(x)
        combs <- expand.grid(
          upstream = 1:(n-1),
          downstream = 2:n
        )
        combs <- combs[combs$downstream > combs$upstream, ]
        ddelta <- x[combs$downstream] - x[combs$upstream]
        idx <- which.max(abs(ddelta))
        tibble(
          MPRO_comparsion = paste0(feature_id[combs$upstream[idx]],"_vs_", feature_id[combs$downstream[idx]]),
          MPRO = ddelta[idx]
        )
    })
  return(MPRO)
}

methods <- c("gene","ALE","TUTR")
for (method in methods){
if (method=="gene") {
  message("Differential APA gene analysis")
    counts=count[,c("gene_name", "PAS_ID",sample_id,"UTR_id")]
    colnames(counts)[1] <- "gene_id"
    colnames(counts)[2] <- "feature_id"
    counts$strand <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[3])
    counts$start <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[2])
    counts <- counts %>%
      mutate(group_id = ifelse(is.na(UTR_id) | UTR_id == ".", feature_id, UTR_id)) %>%
      group_by(gene_id, group_id, strand) %>%
      summarise(
        across(all_of(sample_id), sum),
        start = first(start),  
        feature_id = first(group_id),
        .groups = "drop"
      ) %>%select(-group_id)
    counts <- counts %>% filter(rowSums(across(all_of(sample_id))) > 10)
    merged_data <- counts %>% 
      select(gene_id,feature_id) %>%
      group_by(gene_id) %>%
      summarize(feature_id = str_c(feature_id, collapse = ";"))
    counts2 <- as.data.frame(counts)
    counts <- as.data.frame(counts[,c("feature_id",sample_id,"gene_id")])
    counts <- counts %>%
      distinct(feature_id, .keep_all = TRUE)
  }  
  else if(method=="ALE"){
   message("Differential ALE analysis")
  counts <- count[count$Feature=="3UTR",]
  counts=counts[,c("gene_name", "UTR_id", sample_id)]
  counts <- counts %>% filter(rowSums(across(all_of(sample_id))) > 10)
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

}else if(method=="TUTR"){
  message("Differential TUTR analysis")
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
d <- dmFit(d, design = design_full, verbose = 1)
testgroup <- paste0("Group",group2_sample)
d <- dmTest(d, coef=testgroup)
res <- DRIMSeq::results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]
group1 <- sample[sample$Group==group1_sample,"sample_id"]
group2 <- sample[sample$Group==group2_sample,"sample_id"]

if(method=="TUTR"){
counts$start <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[2])
counts$strand <- sapply(strsplit(as.character(counts$feature_id), ":"), function(x) x[3])
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
MPRO <- calculate_MPRO(counts)
counts$start <- as.numeric(counts$start)
counts_rel <- counts %>%
  group_by(gene_id) %>%
  mutate(
    relative_pos = case_when(
      strand[1] == "+" ~
        (start - min(start)) / (max(start) - min(start)),
      
      strand[1] == "-" ~
        (max(start) - start) / (max(start) - min(start))
    )
  ) %>%
  ungroup()
count_data <- counts_rel %>%
  select(-feature_id, -start, -strand) %>%
  pivot_longer(
    cols = all_of(sample_id),
    names_to = "sample_id",
    values_to = "Count"
  )
count_data <- count_data %>% group_by(gene_id, sample_id) %>% mutate(PAU = Count / sum(Count)) %>% ungroup() 
count_data$gPAU <- count_data$PAU * count_data$relative_pos
WARM_matrix <- count_data %>%
  group_by(gene_id, sample_id) %>%
  summarize(WARM = calculate_WARM(gPAU), .groups = 'drop') %>%
  pivot_wider(names_from = sample_id, values_from = WARM)
WARM_matrix <- as.data.frame(WARM_matrix)
WARM_matrix$WARM1 <-  rowMeans(WARM_matrix[group1],na.rm = TRUE)
WARM_matrix$WARM2 <-  rowMeans(WARM_matrix[group2],na.rm = TRUE)
WARM_matrix$dWARM <- WARM_matrix$WARM1 -WARM_matrix$WARM2
}else if (method=="gene"){
  counts2 <- counts2 %>%
    arrange(gene_id, 
            case_when(
              strand == '-' ~ desc(as.integer(start)), 
              TRUE ~ as.integer(start)  
            )
    )
  unique_sites <- counts2 %>% select(gene_id, feature_id) %>% distinct()
  gene_site_counts <- unique_sites %>% group_by(gene_id) %>%  summarize(Site_Count = n())
  filtered_genes <- gene_site_counts %>% filter(Site_Count > 1) %>% pull(gene_id)
  counts2 <- counts2 %>% filter(gene_id %in% filtered_genes)
  MPRO <- calculate_MPRO(counts2)
  count_data <- counts2 %>%  select(-feature_id, -start, -strand) %>%  pivot_longer(-gene_id, names_to = 'sample_id', values_to = 'Count')
  count_data <- count_data %>% group_by(gene_id, sample_id) %>% mutate(PAU = Count / sum(Count)) %>% ungroup()
  WARM_matrix <- count_data %>%
    group_by(gene_id, sample_id) %>%
    summarize(WARM = calculate_gDPAU(PAU), .groups = 'drop') %>%
    pivot_wider(names_from = sample_id, values_from = WARM)
  WARM_matrix <- as.data.frame(WARM_matrix)
  WARM_matrix$WARM1 <-  rowMeans(WARM_matrix[group1],na.rm = TRUE)
  WARM_matrix$WARM2 <-  rowMeans(WARM_matrix[group2],na.rm = TRUE)
  WARM_matrix$dWARM <- WARM_matrix$WARM1 -WARM_matrix$WARM2
}
else if (method=="ALE"){
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
  MPRO <- calculate_MPRO(counts)
  count_data <- counts %>%  select(-feature_id, -utr) %>%  pivot_longer(-gene_id, names_to = 'sample_id', values_to = 'Count')
  count_data <- count_data %>% group_by(gene_id, sample_id) %>% mutate(PAU = Count / sum(Count)) %>% ungroup()
  
  WARM_matrix <- count_data %>%
    group_by(gene_id, sample_id) %>%
    summarize(WARM = calculate_gDPAU(PAU), .groups = 'drop') %>%
    pivot_wider(names_from = sample_id, values_from = WARM)
  WARM_matrix <- as.data.frame(WARM_matrix)
  WARM_matrix$WARM1 <-  rowMeans(WARM_matrix[group1],na.rm = TRUE)
  WARM_matrix$WARM2 <-  rowMeans(WARM_matrix[group2],na.rm = TRUE)
  WARM_matrix$dWARM <- WARM_matrix$WARM1 -WARM_matrix$WARM2
}else {
  message("not available")
}
res <- res %>% left_join(WARM_matrix, by ="gene_id")
res <- res %>% left_join(MPRO,by ="gene_id")
res <- res %>% left_join(merged_data, by ="gene_id")
res$Note <- ifelse(res$MPRO > 0, "distal", "proximal")
res$sig <- ifelse(abs(res$MPRO) >= 0.2 & res$adj_pvalue <= 0.05 , "TRUE", "FALSE")
write.table(res, paste0(group1_sample, '_', group2_sample, '.', method, ".diff.DRIMSeq.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}
