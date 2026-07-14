---
layout: default
title: Step_4_Diff
nav_order: 7
---

# **Step 4. Differential analysis**

![](images/PAS.png)

## PAS usage index

We quantified PAs using the metric “poly(A) site usage” (PAU). The PAU for a PA is the ratio of its read count to the sum of read counts of all detected PAs from its gene, which ranges from 0 to 1. We used the weighted average relative mode (WARM) value to quantify the proximal-distal relative PAU of each gene or 3'-UTR as [we previous described](https://genome.cshlp.org/content/33/10/1774) with some minor modifications. For a gene with n (n≥2) PAs: (1) if all PAs were annotated within same 3′ UTR, the most proximal and most distal PAs were assigned relative positions of 0 and 1, respectively, and the relative positions of intermediate PAs were linearly interpolated according to their genomic coordinates. (2) if PAs located in different 3′ UTRs and/or non-3′ UTRs, PAs were ranked from 0 to (n-1) according to their genomic position. The relative position of the PA rank *i* was defined as *i*/(n-1). WARM was calculated as the weighted average of relative PA positions using PAU values as weights. WARM value ranges from 0 to 1, where higher values indicate a greater likelihood of distal PA usage within a gene or a 3'-UTR, and *vice versa*.

The [maximum difference in proportion change (MPRO) value](https://genome.cshlp.org/content/33/10/1774) was used to measure the largest shift in relative PA usage between two groups. The PAU difference (δ) of each PA between the base and alternative (alt) groups was first calculated. For a PA pair (i, j), where PA i is located downstream of PA j, the difference in proportion change was defined as dδ = δ~i~ – δ~j~. MPRO was defined as the dδ value with the largest absolute magnitude across all possible PA pairs. A positive or negative MPRO value indicates increased distal or proximal PA usage, respectively, in the base group relative to the alt group.  

## 4.1. With replicates

``` shell
lrapa diff -i test.PAS.count.txt -s sample.txt -c case -n control
```

This module performs differential APA usage analyses between exactly two conditions with two or more replicates based the R package [DRIMSeq](http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html). This is done by testing if the ratio of APA changes between conditions.

If a gene's absolute mean difference of MPRO is \>= 0.2 and adjusted *P* value is \<= 0.05 between two groups, the gene's APA change between two groups will be deemed as significant.

**Parameters**

- `-i, --input`: count file of PAS sites (required).
- `-s, --sample`: sample information (required).
- `-c, --group1`: group1 for differential analysis (required).
- `-n, --group2`: group2 for differential analysis (required).

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

If a gene's absolute mean difference of MPRO is \>= 0.4 and adjusted *P* value is \<= 0.05 between two groups, APA change between two groups will be deemed as significant.

``` shell
lrapa norepdiff -i test.PAS.count.txt -c sample1 -n sample2
```

**Parameters**

- `-i, --input`: count file of PAS sites (required).

- `-c, --group1`: sample 1 for differential analysis (required).

- `-n, --group2`: sample 2 for differential analysis (required).

  **Note:** The name of sample must same as in the count file.

The diff function categorizes genes into three classes based on differential usage patterns: (1) Differential APA genes: Genes in which multiple PAs are differentially utilized, non-3'-UTR PAs are also included for statistical analysis; (2) Differential TUTR-APA genes: Genes where multiple PAs within the same 3'-UTR are differentially utilized; (3) Differential ALS-APA genes: Genes in which multiple 3'-UTRs are differentially utilized, combining multiple PAs within the 3'-UTR. The outputs generated three files for these differential usage patterns.
