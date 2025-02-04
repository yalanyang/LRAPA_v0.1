---
layout: default
title: Step_3_Count
nav_order: 6
---

# **Step 3. PAS count**

## 3.1 Bulk long-read RNA-seq

``` shell
lrapa count -b N1.bam,N2.bam,C1.bam,C2.bam -p test.PAS.bed -o test.PAS.count.txt
```

**Parameters**

-   `-b, --bam_files`: Comma-separated list of input BAM files (required).
-   `-p, --pas`: PAS reference TXT file with BED columns followed by annotations (required).
-   `-o, --output`: Output text file for the count matrix (required).

The output count format:

|                                  | C1  | C2  | C3  | N1  | N2  | N3  | Feature | gene_name |             UTR_id              |       PAS_ID        |
|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
| chr11_35229465_35229470\_+\_3UTR | 16  | 30  | 22  | 12  | 10  |  8  |  3UTR   |   CDC44   | CD44:chr11:35229288:35232402:u2 | chr11_35229470_3UTR |
| chr11_35230014_35230028\_+\_3UTR | 69  | 50  | 54  | 34  | 32  | 22  |  3UTR   |   CDC44   | CD44:chr11:35229288:35232402:u2 | chr11_35230028_3UTR |
| chr11_35232387_35232403\_+\_3UTR | 140 | 160 | 155 | 310 | 440 | 550 |  3UTR   |   CDC44   | CD44:chr11:35229288:35232402:u2 | chr11_35232403_3UTR |

## 3.2 Single-cell long-read RNA-seq

``` shell
lrapa scCount -b scLR.bam -p test.PAS.bed -o test -w barcode.txt
```

**Parameters**

-   `-i, --input`: Input scRNA-seq BAM file (required).
-   `-p, --pas`: PAS reference TXT file with BED columns followed by annotations (required).
-   `-o, --output`: Output folder for the count matrix (required).
-   `-w, --barcode`: Barcode file.