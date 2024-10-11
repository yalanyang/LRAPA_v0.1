---
layout: default
title: Step_1_Filtering
nav_order: 4
---


# **Step 1. Filter long-read reads with polyA signals**

> To identify Iso-Seq reads that capture cleavage and polyadenylation events, we searched for reads that contained stretches of adenosines (i.e., polyA tails).
>
> The poly(A) cleavage site on the genome is considered to be right after the 3′-most position of the alignment of long reads with the genome.
>
> It contains a soft-clipped sequence mostly composed of adenines (or thymines for reverse strand). PolyA stretches needed to be located immediately after the cleavage site (i.e., starts at the base within the read that does not map to the hg38 reference genome, which is otherwise known as the portion of the read that is "softclipped"). We assessed if the softclipped portion of every read contained a stretch of adenosines. We retained the reads if their softclipped segments were \< 20 nucleotides in length and were composed of 95% adenosines. Moreover, if the length of the softclipped segment of a read was ≥ 20 nucleotides, we assessed if the first 20 nucleotides of the softclipped segment was composed of 80% adenosines and if the following 20 nucleotides of the softclipped segment was composed of 95% adenosines and retained these reads as containing a stretch of adenosines.
>
> The flanking region does not contain a stretch of six adenines. Reads with stretches of adenosines were filtered for internal priming or mispriming using an approach similar to what has been described previously. In brief, the genomic sequence −10 to +10 nt surrounding the cleavage site was examined(`-f, --flank_size)`. If the sequence has six continuous As, it is considered as internal priming.
>
> The adjacent region contains a known PAS hexamer. We require the long- read has one of the 12 PAS hexamers (AAUAAA or 11 variants) in −40 to −1 nt region of the cleavage site (`-a, --adjacent_size`).

The 12 PAS hexamers used here are "AATAAA", "TTTAAA", "AAGAAA", "AACAAA", "TATAAA", "AATGAA", "ATTAAA", "AGTAAA", "AATATA", "CATAAA", "ACTAAA", "GATAAA", which is adopted from [our previous study.](https://genome.cshlp.org/content/33/10/1774.full)

``` shell
lrapa filter -i test.flnc.filter.bam -r $reference/GRCh38.primary_assembly.genome.fa -o test.HQ.qname.txt
```

**Parameters**

-   `-i, --input`: Input mapped BAM file (required).
-   `-o, --hq`: Output BAM file (default: "HQ.bam").
-   `-r, --reference`: Reference genome FASTA file (required).
-   `-f, --flank_size`: Flank size on either side of the mapped end position for internal priming (default: 10).
-   `-c, --cores`: Number of cores to use for parallel processing (default: 4).
-   `-a, --adjacent_size`: Flank size to search for PAS signal (default: 40).
-   `-b, --batch_size`: Number of reads to process in each batch (default: 10000).

The script will output a txt file (test.HQ.qname.txt) contained the qname of the reads with ployA signals. Then we extract these reads from the bam file with samtools. The extract reads will be used for PAS identification and annotation analysis in step 2.

``` shell
samtools view -H test.fl_flipped.unique.bam > test.header
samtools view -N test.HQ.qname.txt -o test.HQ.sam test.fl_flipped.unique.bam
cat test.header test.HQ.sam > test.HQ.header.sam 
samtools view -bS test.HQ.header.sam > test.HQ.bam
```

**Note:** Make sure the version of samtools is 1.12 or greater, which accepts option `-N` . If you find that this script needs a lot of memory or very slow you may want to split the input bam file by chromosome (samtools view) and run these separately. We do intend to improve this.
