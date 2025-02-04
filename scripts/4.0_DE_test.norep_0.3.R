# Load necessary libraries

library(optparse)

option_list <- list(
  make_option(c("-i", "--count_files"), type = "character", default = NULL, 
              help = "Comma-separated list of input count files", metavar = "character"),
  make_option(c("-c", "--group1"), type = "character", default = NULL, 
              help = "Comma-separated list of input count files", metavar = "character"),
  make_option(c("-n", "--group2"), type = "character", default = NULL, 
              help = "Comma-separated list of input count files", metavar = "character")
)

##Usage:  Rscript 4.0_DE_test.norep_0.2.R -i Encode.brain.bulk.count.txt -c AD -n Health

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

library(dplyr)
count <-opt$count_files

PAS <-  read.table(count, sep = "\t", header = T, stringsAsFactors = FALSE)

group1 <- opt$group1
group2 <- opt$group2

Count1_total <- aggregate(PAS[,group1],by=list(type=PAS$gene_name),sum)
colnames(Count1_total) <- c("gene_name","gene_count1")
Count2_total <- aggregate(PAS[,group2],by=list(type=PAS$gene_name),sum)
colnames(Count2_total) <- c("gene_name","gene_count2")


PAS = PAS %>% left_join(Count1_total, by="gene_name")
PAS = PAS %>% left_join(Count2_total, by="gene_name")
PAS$PAU1 <- round(PAS[,group1]/PAS$gene_count1,3)
PAS$PAU2 <- round(PAS[,group2]/PAS$gene_count2,3)
PAS$dPAU <- PAS$PAU1-PAS$PAU2

PAS <- PAS %>% dplyr::filter(((PAU1>=0.05 | PAU2>=0.05) & (PAU1 <=0.95 | PAU2 <= 0.95)) & (PAS[,group1] >=10 | PAS[,group2]>=10))


# Calculate gDPAU for each gene and strand
calculate_gDPAU <- function(df) {
  n <- length(df)
  p <- df
  gDPAU <- sum((seq_len(n) - 1) * p) / (n - 1)
  return(gDPAU)
}

betweenPAS <- function(data) {
  # Perform chi-squared tests for each gene
  genes <- unique(data$gene_name)
  results <- data.frame(
    Gene = character(),
    Comparsions = character(),
    UTR_id = character(),
    p_value = numeric(),
    group1_count =character(),
    group2_count = character(),
    group1_PAU =  character(),
    group2_PAU =  character(),
    PDUI <- numeric(),
    Note1 =  character(),
    Note2 =  character(),
    stringsAsFactors = FALSE
  )
  
  for (gene in genes) {
    gene_data <- subset(data, gene_name == gene)
    # Get all pairs of rows for the gene
    num_sites <- nrow(gene_data)
    if (num_sites > 1) {
      
      for (i in 1:(num_sites-1)) {
        for (j in (i+1):num_sites) {
          pair <- gene_data[c(i, j), ]
          contingency_table <- pair[, c(group1, group2)]
          chisq_test <- chisq.test(contingency_table)
          pair_description <- paste0(gene_data$PAS_ID[i], " vs ", gene_data$PAS_ID[j])
          UTR <- paste0(gene_data$UTR_id[i], " vs ", gene_data$UTR_id[j])
          group1_count <- paste0(gene_data[,group1][i], ",", gene_data[,group1][j])
          group2_count <- paste0(gene_data[,group2][i], ",", gene_data[,group2][j])
          group1_PAU <- paste0(gene_data$PAU1[i],",",gene_data$PAU1[j])
          group2_PAU <- paste0(gene_data$PAU2[i],",",gene_data$PAU2[j])
          strand <- strsplit(gene_data$PAS_ID[i],"[:]")[[1]][3]
          start_group1 <- as.numeric(strsplit(gene_data$PAS_ID[i],"[:]")[[1]][2])
          start_group2  <- as.numeric(strsplit(gene_data$PAS_ID[j],"[:]")[[1]][2])
          
          DPAU <- ifelse ((strand == "+" & start_group1 > start_group2) | (strand == "-" & start_group1 < start_group2),  gene_data$PAU1[i] - gene_data$PAU2[i],
                          gene_data$PAU1[j] - gene_data$PAU2[j])
          
          
          
          if(gene_data$Feature[i]=="3UTR" & gene_data$Feature[j]=="3UTR")  { 
            
            if (gene_data$UTR_id[i]==gene_data$UTR_id[j]){
              Note1 <- "within 3UTR"
              Note2 <- ifelse(DPAU > 0, "lengthen", "shorten")
            }else {
              Note1 <- "between 3UTR"
              Note2 <- ifelse(DPAU > 0, "Distal UTRs", "Proximal UTRs")
            }
            
          }else {
            Note1 <- paste0(gene_data$Feature[i], ",", gene_data$Feature[j])
            Note2 <- "NA"
          }
          
          results <- rbind(results, data.frame(
            Gene = gene_data$gene_name[i],
            Comparsions = pair_description,
            UTR_id = UTR,
            p_value = chisq_test$p.value,
            group1_count = group1_count,
            group2_count = group2_count,
            group1_PAU = group1_PAU,
            group2_PAU = group2_PAU,
            DPAU = DPAU,
            Note1 =   Note1,
            Note2 =   Note2,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  results <- results[order(results$p_value), ]
  results$adj_pvalue <- p.adjust(results$p_value, method = "BH")
  results$sig <- ifelse(abs(results$DPAU) >= 0.1 & results$adj_pvalue <= 0.05 , "TRUE", "FALSE")
  return(results)
}

betweenPAS_gene <- function(data) {
# Perform chi-squared tests for each gene
genes <- unique(data$gene_name)
results <- list()
data$start <- sapply(strsplit(as.character(data$PAS_ID), ":"), function(x) x[2])
data$strand <- sapply(strsplit(as.character(data$PAS_ID), ":"), function(x) x[3])
for (i in seq_along(genes)) {
  count <- subset(data, gene_name == genes[i])
  # Get all pairs of rows for the gene
  num_sites <- nrow(count)
  if (num_sites >= 2) {
    
    contingency_table <- as.matrix(count[, c(group1, group2)])
    chisq_test <- chisq.test(contingency_table)
    pas_list <- paste(count$PAS_ID, collapse = ",")
    count1_list <- paste(count[,group1], collapse = ",")
    count2_list <- paste(count[,group2], collapse = ",")
    group1_PAU <-  paste(count$PAU1, collapse = ",")
    group2_PAU <-  paste(count$PAU2, collapse = ",")
    # PDUI <- ifelse (strand == "+", count[count$start== max(count$start),]$PAU1 - count[count$start==max(count$start),]$PAU2, 
    #                 count[count$start== min(count$start),]$PAU1 - count[count$start==min(count$start),]$PAU2)
    strand <- count$strand[1]
    count <- count %>%
      arrange(
        case_when(
          strand == '-' ~ desc(as.integer(start)),  # Convert start to integer if needed
          TRUE ~ as.integer(start)  # Convert start to integer if needed
        )
      )
    gDPAU1 <- calculate_gDPAU(count$PAU1)
    gDPAU2 <- calculate_gDPAU(count$PAU2)
    
    dgDPAU <- gDPAU1 -gDPAU2
    Note <- ifelse (dgDPAU > 0, "distal", "proximal")
    
    result <- data.frame(
      gene_name = unique(count$gene_name),
      p_value = chisq_test$p.value,
      pas_list = pas_list,
      count1_list = count1_list,
      count2_list = count2_list,
      group1_PAU =   group1_PAU,
      group2_PAU =   group2_PAU,
      gDPAU1 = gDPAU1,
      gDPAU2 = gDPAU2,
      dgDPAU = dgDPAU,
      Note = Note,
      stringsAsFactors = FALSE
    )
    results[[i]] <- result
  }
}
results_df <- do.call(rbind, results)
results_df <- results_df[order(results_df$p_value), ]
results_df$adj_pvalue <- p.adjust(results_df$p_value, method = "BH")
results_df$sig <- ifelse(abs(results_df$dgDPAU) >= 0.1 & results_df$adj_pvalue <= 0.05 , "TRUE", "FALSE")
return(results_df)
}

within_3UTR <- function(data) {
  results <- list()
  data$start <- sapply(strsplit(as.character(data$PAS_ID), ":"), function(x) x[2])
  data$strand <- sapply(strsplit(as.character(data$PAS_ID), ":"), function(x) x[3])
  utrs <- unique(data$UTR_id)
  for ( i in seq_along(utrs)){
    count <- data[data$UTR_id==utrs[i],]
    if (nrow(count) > 1) {
      contingency_table <- as.matrix(count[, c(group1, group2)])
      chisq_test <- chisq.test(contingency_table)
      pas_list <- paste(count$PAS_ID, collapse = ",")
      count1_list <- paste(count[,group1], collapse = ",")
      count2_list <- paste(count[,group2], collapse = ",")
      group1_PAU <-  paste(count$PAU1, collapse = ",")
      group2_PAU <-  paste(count$PAU2, collapse = ",")
      strand <- count$strand[1]
     # PDUI <- ifelse (strand == "+", count[count$start== max(count$start),]$PAU1 - count[count$start==max(count$start),]$PAU2, 
     #                 count[count$start== min(count$start),]$PAU1 - count[count$start==min(count$start),]$PAU2)
      
      count <- count %>%
        arrange(
          case_when(
            strand == '-' ~ desc(as.integer(start)),  # Convert start to integer if needed
            TRUE ~ as.integer(start)  # Convert start to integer if needed
          )
        )
      
      gDPAU1 <- calculate_gDPAU(count$PAU1)
      gDPAU2 <- calculate_gDPAU(count$PAU2)
      
      dgDPAU <- gDPAU1 -gDPAU2
      Note <- ifelse (dgDPAU > 0, "lengthen", "shorten")
      
      result <- data.frame(
        utr_id = unique(count$UTR_id),
        gene_name = unique(count$gene_name),
        p_value = chisq_test$p.value,
        pas_list = pas_list,
        count1_list = count1_list,
        count2_list = count2_list,
        group1_PAU =   group1_PAU,
        group2_PAU =   group2_PAU,
        gDPAU1 = gDPAU1,
        gDPAU2 = gDPAU2,
        dgDPAU = dgDPAU,
        Note = Note,
        stringsAsFactors = FALSE
      )
      results[[i]] <- result
    }
  }
  results_df <- do.call(rbind, results)
  # Sort results by p-value in ascending order
  results_df <- results_df[order(results_df$p_value), ]
  #results_df$q_value <- p.adjust(results_df$p_value, method = "bonferroni")
  results_df$adj_pvalue <- p.adjust(results_df$p_value, method = "BH")
  results_df$sig <- ifelse(abs(results_df$dgDPAU) >= 0.1 & results_df$adj_pvalue <= 0.05 , "TRUE", "FALSE")
  return(results_df)
}


ALE <- function(data) {
  results <- list()
  data$strand <- sapply(strsplit(as.character(data$PAS_ID), ":"), function(x) x[3])
  Count1_total <- aggregate(data[,group1],by=list(type=data$UTR_id),sum)
  PAU1_total <-  aggregate(data$PAU1,by=list(type=data$UTR_id),sum)
  colnames(Count1_total) <- c("UTR_id","utr_count1")
  colnames(PAU1_total) <- c("UTR_id","PAU1")
  Count2_total <- aggregate(data[,group2],by=list(type=data$UTR_id),sum)
  PAU2_total <-  aggregate(data$PAU2,by=list(type=data$UTR_id),sum)
  colnames(Count2_total) <- c("UTR_id","utr_count2")
  colnames(PAU2_total) <- c("UTR_id","PAU2")
  data = unique(data %>% select("UTR_id","gene_name","strand"))
  data = data %>% left_join(Count1_total, by="UTR_id") %>% 
                  left_join(Count2_total, by="UTR_id") %>% 
                  left_join(PAU1_total, by="UTR_id") %>% 
                  left_join(PAU2_total, by="UTR_id")
  data$start <- sapply(strsplit(as.character(data$UTR_id), ":"), function(x) x[4])
  genes <- unique(data$gene_name)
  for ( i in seq_along(genes)){
    count <- data[data$gene_name==genes[i],]
    count <- count %>% filter((utr_count1+utr_count2)>=10)
    if (nrow(count) > 1) {
      contingency_table <- as.matrix(count[, c("utr_count1","utr_count2")])
      chi2_test <- chisq.test(contingency_table)
      utr_list <- paste(count$UTR_id, collapse = ",")
      count_group1 <- paste(count$utr_count1, collapse = ",")
      count_group2 <- paste(count$utr_count2, collapse = ",")
      PAU_group1 <- paste(count$PAU1, collapse = ",")
      PAU_group2 <- paste(count$PAU2, collapse = ",")
      strand <- count$strand[1]
      count <- count %>%
        arrange(
          case_when(
            strand == '-' ~ desc(as.integer(start)),  # Convert start to integer if needed
            TRUE ~ as.integer(start)  # Convert start to integer if needed
          )
        )
      
      gDPAU1 <- calculate_gDPAU(count$PAU1)
      gDPAU2 <- calculate_gDPAU(count$PAU2)
      
      dgDPAU <- gDPAU1 -gDPAU2
      
      #PDUI <- ifelse (strand == "+", count[count$start== max(count$start),]$PAU1 - count[count$start==max(count$start),]$PAU2, 
       #               count[count$start== min(count$start),]$PAU1 - count[count$start==min(count$start),]$PAU2)
      Note <- ifelse (dgDPAU > 0, "distal", "proximal")
 
      result <- data.frame(
        gene_name = unique(count$gene_name),
        p_value = chi2_test$p.value,
        utr_list = utr_list,
        count_group1 = count_group1,
        count_group2 = count_group2,
        PAU_group1 =  PAU_group1,
        PAU_group2 =  PAU_group2,
        gDPAU1 = gDPAU1,
        gDPAU2 = gDPAU2,
        dgDPAU = dgDPAU,
        Note = Note,
        strand = strand,
        stringsAsFactors = FALSE
      )
      results[[i]] <- result
    }
  }
  results_df <- do.call(rbind, results)
  results_df <- results_df[order(results_df$p_value), ]
  #results_df$q_value <- p.adjust(results_df$p_value, method = "bonferroni")
  results_df$adj_pvalue <- p.adjust(results_df$p_value, method = "BH")
  results_df$sig <- ifelse(abs(results_df$dgDPAU) >= 0.1 & results_df$adj_pvalue <= 0.05 , "TRUE", "FALSE")
  return(results_df)
}


# Perform chi-squared test and adjust p-values
between_PAS_gene <- betweenPAS_gene(PAS)
write.table(between_PAS_gene, paste0(group1, '_', group2, ".diff.betweenPAS_gene_level.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

#between_PAS <- betweenPAS(PAS)
#write.table(between_PAS, paste0(group1, '_', group2, ".diff.betweenPAS.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

PAS <- PAS[PAS$Feature=="3UTR",]
within_3UTR <- within_3UTR(PAS)
between_3UTR <- ALE(PAS)
write.table(within_3UTR, paste0(group1, '_', group2, ".PAS.diff.TUTR.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(between_3UTR, paste0(group1, '_', group2, ".PAS.diff.ALS.txt"), sep = "\t", row.names = FALSE, quote = FALSE)   
cat("Analysis complete", "\n")
