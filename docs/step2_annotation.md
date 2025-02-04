---
layout: default
title: Step_2_Annotation
nav_order: 5
---

# **Step 2. PAS identification and annotation**

> After obtaining long reads containing poly(A) signals, we identify the poly(A) cleavage site for each read, defined as the last mapped base of the read. Since the cleavage can be imprecise, resulting in mRNAs with variable ends, we refer to the cleavage site as a location where mRNA cleavage takes place, and poly(A) site as a region containing cleavage site(s). Here, due to the inherent heterogeneity of polyadenylation cleavage, we iteratively clustered poly(A) cleavage sites that are within 24 nucleotides of each other.
>
> If a poly(A) site contained a cleavage site annotated in the reference GTF file, the annotated cleavage site is used as the representative poly(A) site. If no annotated site within the poly(A) site, we select the cleavage site with an adenine ("A") base and the highest read count as the representative poly(A) site. If none of the poly(A) cleavage sites are an "A" base, we choose the cleavage site with the highest read coverage as the representative poly(A) site.
>
> **reference**: Bin Tian, Jun Hu, Haibo Zhang, Carol S. Lutz, A large-scale analysis of mRNA polyadenylation of human and mouse genes, *Nucleic Acids Research*, Volume 33, Issue 1, 1 January 2005, Pages 201--212. <https://doi.org/10.1093/nar/gki158>

``` r
lrapa anno -i test.HQ.bam -r genome.fa -g refGene.gtf -u refGene.3utr.bed -o test.PAS.bed
```

**Note:** The hg38.refGene.3utr_merge.bed was generated using the get3UTR_0.1.R script. if you find that this script needs a lot of memory or very slow you may want to split the input bam file by chromosome and run these separately. We do intend to improve this.

**Parameters**

-   `-i, --input`: Input BAM file (required).
-   `-r, --reference`: Reference genome FASTA file (required).
-   `-g, --gtf`: Reference annotation gtf file (required).
-   `-f, --flank_size`: Flank size to search for PAS signal (default: 40).
-   `-u, --utr`: Reference annotation 3'-UTR file (required).
-   `-o, -–output`: Output PAS bed file with header (required).