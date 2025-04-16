---
layout: default
title: Step_1_Filtering
nav_order: 4
---

# **Step 1. Filter long-read reads with polyA signals**

This step is used to obtain long-reads capturing cleavage and polyadenylation events.

> To obtain reads capturing cleavage and polyadenylation events, the 3′-most position of mapped base of the read was tentatively considered as the putative cleavage site. The 3′-end of the read was required to contain a soft-clipped sequence, mostly composed of adenines (or thymines for the reverse strand to ensure the read contain a polyA tail. Reads were retained if the soft-clipped sequence was \<20 nucleotides and ≥95% adenosines, or if ≥20 nucleotides, the first 20 contained ≥80% adenosines, and the next 20 had ≥ 95%. To exclude internal priming, the −10 to +10 nt region around the putative cleavage site could not contain a stretch of six consecutive adenines. Additionally, the region −40 to −1 nt upstream of the putative cleavage site had to include one of 12 canonical PA hexamers (AAUAAA, AUUAAA, UUUAAA, AAGAAA, AACAAA, UAUAAA, AAUGAA, AGUAAA, AAUAUA, CAUAAA, ACUAAA, GAUAAA).

``` shell
lrapa filter -i test.flnc.filter.bam -r ref.genome.fa -o test.HQ.qname.txt
```

**Parameters**

-   `-i, --input`: Input mapped BAM file (required).
-   `-o, --hq`: Output BAM file (default: "HQ.bam").
-   `-r, --reference`: Reference genome FASTA file (required).
-   `-f, --flank_size`: Flank size on either side of the mapped end position for internal priming (default: 10).
-   `-c, --cores`: Number of cores to use for parallel processing (default: 4).
-   `-a, --adjacent_size`: Flank size to search for PAS signal (default: 40).
-   `-b, --batch_size`: Number of reads to process in each batch (default: 10000).

The script will output a txt file (test.HQ.qname.txt) contained the qname of the reads with ployA signals. Then we extract these reads from the bam file with samtools. The extract reads will be used for PA identification and annotation analysis in step 2.

``` shell
samtools view -H test.fl_flipped.unique.bam > test.header
samtools view -N test.HQ.qname.txt -o test.HQ.sam test.fl_flipped.unique.bam
cat test.header test.HQ.sam > test.HQ.header.sam 
samtools view -bS test.HQ.header.sam > test.HQ.bam
```

**Note:** Make sure the version of samtools is 1.12 or greater, which accepts option `-N` . If you find that this script needs a lot of memory or very slow you may want to split the input bam file by chromosome (samtools view) and run these separately. We do intend to improve this.
