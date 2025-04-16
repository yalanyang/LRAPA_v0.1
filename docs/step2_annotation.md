---
layout: default
title: Step_2_Annotation
nav_order: 5
---

# **Step 2. PA identification and annotation**

After filtering (step 1) to obtain the full-length reads, the BAM files from all samples should be merged using Samtools (samtools merge). After merging, the next step is polyadenylation (PA) site identification and annotation.

``` r
lrapa anno -i test.HQ.bam -r ref.genome.fa -g refGene.gtf -u refGene.3utr.bed -o test.PAS.bed
```
**Parameters**

-   `-i, --input`: Input BAM file (required).
-   `-r, --reference`: Reference genome FASTA file (required).
-   `-g, --gtf`: Reference annotation gtf file (required).
-   `-f, --flank_size`: Flank size to search for poly(A) signal (default: 40).
-   `-u, --utr`: Reference annotation 3'-UTR file (required).
-   `-o, -–output`: Output PAS bed file with header (required).

**Note:** if you find that this script needs a lot of memory or very slow you may want to split the input bam file by chromosome and run these separately. We do intend to improve this.

The 3'-UTR reference (refGene.3utr.bed) is generated using the get3UTR_0.1.R script in the scripts file. 

