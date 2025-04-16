---
layout: default
title: Step_2_Annotation
nav_order: 5
---

# **Step 2. PA identification and annotation**

After filtering (step 1) to obtain the full-length reads, the BAM files from all samples should be merged using Samtools (samtools merge). After merging, the next step is polyadenylation (PA) site identification and annotation.

> We define a cleavage site as the specific base location where mRNA cleavage occurs and cluster cleavage sites within 24 nucleotides of each other into a PA, which is a region that may contain one or more cleavage sites. Within each PA, the representative cleavage site was determined as follows: if an annotated cleavage site was present in the reference, it is selected. If no annotated site was available, the cleavage site with an adenine base and the highest read count was chosen. If none of the cleavage sites contained an adenine, the site with the highest read coverage is selected as the representative cleavage site. PAs are annotated according to RefGene. For PAs falling into multiple annotated regions, we set priority as 3′-UTR \> stop codon \> 5′-UTR \> start codon \> exon \> intron \> intergenic.

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

The 3'-UTR reference (refGene.3utr.bed) is generated using the get3UTR_0.1.R script in the scripts file, as shown in the [**Other**](https://yalanyang.github.io/LRAPA_v0.1/others.html) section.
