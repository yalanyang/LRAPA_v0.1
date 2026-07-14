---
layout: default
title: Step_4_Diff
nav_order: 7
---

# **Step 4. Differential analysis**

![](images/PAS.png)

## PAS usage index


We quantified PAs using the metric “polyA site usage” (PAU). The PAU for a PA is the ratio of its read count to the sum of read counts of all detected PAs from its gene, which ranges from 0 to 1. 

For genes with two polyA sites, we used DPAU to quantify the percentage of DPAU for each gene. DPAU ranges from 0 to 1.

![](images/DPAU.png)

For genes or 3’-UTRs with two PAs, we used the distal polyA site usage index (DPAU) to quantify the PAU usage of each gene or 3'-UTR, which is the ratio of reads at distal PAs to the sum of reads at distal and proximal PAs. 

For genes or 3'-UTRs with more than two PAs, generalized DPAU ([gDPAU](https://www.pnas.org/doi/10.1073/pnas.2113504119) was used to quantify the trend of distal PA usage for each gene or 3’-UTR. It is a location index weighted PAU of each PA. When n = 2, DPAU = gDPAU. DPAU and gDPAU ranges from 0 to 1, where higher values indicate a greater likelihood of distal PA usage within a gene or a 3'-UTR, and vice versa. 

![](images/gDUAP.png)

When n = 2, gDPAU = DPAU.

## 4.1. With replicates

``` shell
lrapa diff -i test.PAS.count.txt -s sample.txt -c case -n control
```

This module performs differential APA usage analyses between exactly two conditions with two or more replicates based the R package [DRIMSeq](http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html). This is done by testing if the ratio of APA changes between conditions.

If a gene's absolute mean difference of DPAU/gDPAU is \>0.1 and adjusted *P* value is \< 0.05 between two groups, the gene's APA change between two groups will be deemed as significant.

**Parameters**

-   `-i, --input`: count file of PAS sites (required).
-   `-s, --sample`: sample information (required).
-   `-c, --group1`: group1 for differential analysis (required).
-   `-n, --group2`: group2 for differential analysis (required).

Example file of sample information

| sample_id,group |
|:---------------:|
|     N1,case     |
|     N2,case     |
|     N3,case     |
|   C1,control    |
|   C2,control    |
|   C3,control    |

## 4.2. No replicates

This module performs differential APA usage analyses between exactly two conditions with no replicate based chisq test. For each test, a n × 2 matrix per gene or 3'-UTR was generated, with the n PA forming rows and read count in the two conditions forming columns.

We used a Benjamini--Hochberge correction for multiple testing.

If a gene's absolute mean difference of DPAU/gDPAU is \>0.2 and adjusted *P* value is \< 0.05 between two groups, APA change between two groups will be deemed as significant.

``` shell
lrapa norepdiff -i test.PAS.count.txt -c sample1 -n sample2
```

**Parameters**

-   `-i, --input`: count file of PAS sites (required).
-   `-c, --group1`: sample 1 for differential analysis (required).
-   `-n, --group2`: sample 2 for differential analysis (required).

    **Note:** The name of sample must same as in the count file.

The diff function categorizes genes into three classes based on differential usage patterns: (1) Differential APA genes: Genes in which multiple PAs are differentially utilized, non-3'-UTR PAs are also included for statistical analysis; (2) Differential TUTR-APA genes: Genes where multiple PAs within the same 3'-UTR are differentially utilized; (3) Differential ALS-APA genes: Genes in which multiple 3'-UTRs are differentially utilized, combining multiple PAs within the 3'-UTR. The outputs generated three files for these differential usage patterns.