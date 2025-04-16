---
layout: default
title: Step_0_Pre-processing
nav_order: 3
---

# **Step 0. Pre-processing of the long-read CCS reads**

## 0.1. bulk long-read RNA-seq

LRAPA accept the PacBio CCS bam file as the input. The test data could be found in the test_data files.
Removal of primers from the CCS reads is performed using [lima](https://isoseq.how).

``` shell
lima test.ccs.bam IsoSeq_v2_primers_12.fasta test.fl.bam --isoseq --peek-guess
isoseq3 refine test.fl.IsoSeqX_bc04_5p--IsoSeqX_3p.bam IsoSeq_v2_primers_12.fasta test.refine.bam
```

**Note**: The IsoSeq_v2_primers_12.fasta file was download from [Pacbio](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/)

Then we convert the full-length bam file into sam format and run the [flip_reads.py](https://github.com/mortazavilab/ENCODE-references/tree/master) script to orient the reads to the correct strand (since lima\--ccs does not do this, recommend by [Encode long-read pipeline](https://www.encodeproject.org/documents/6d583a1d-d692-4511-b13b-c051822d861c/@@download/attachment/ENCODE%20Long%20Read%20RNA-Seq%20Analysis%20Pipeline%20v3.2%20%28Human%29.pdf)).

``` shell
samtools view -h test.refine.bam > test.refine.sam
python flip_reads.py --f test.refine.sam --o test.fl_flipped.sam
samtools view -bS test.fl_flipped.sam > test.fl_flipped.bam
isoseq3 refine test.fl_flipped.bam PB_adapters.fasta test.flnc.bam
```

**Note**: We don't trim poly(A) tails because we aim to assess whether each read contains a genuine poly(A) signal, based on the length and base composition of its poly(A) tail.

We next align the long-reads to the reference genome with Minimap2. The human reference genome could found here: [GRCh38 human reference genome](https://www.gencodegenes.org/human/release_21.html) 

``` shell
bamtools convert -format fastq -in test.flnc.bam -out test.flnc.fastq
minimap2  -ax splice -uf -C5 ref.genome.fa test.flnc.fastq > test.flnc.mapping.sam
samtools view -O BAM -F 2052 -h test.flnc.mapping.sam |  samtools sort -O BAM -@ 7 -o test.flnc.unique.bam -
samtools view -h test.flnc.unique.bam | awk '$10 != "*"' |samtools view -bS - > test.flnc.filter.bam
```

## 0.2. Single-cell long-read RNA-seq

For single-cell long-read data, we use the [isoseq pipeline](https://isoseq.how/umi/cli-workflow.html) to process the CCS reads. The tag function clips UMIs and cell barcodes from the reads and associates them with the reads for later deduplication. We retain the poly(A) tail in the isoseq refine function for the downstream APA analysis.

``` shell
lima --isoseq --dump-clips --peek-guess -j 24 test.sc.hifi_reads.bam 10x_Chromium_3p_primers.fasta test.fl.bam
isoseq3 tag test.fl.5p--3p.bam test.5p--3p.tagged.bam --design T-12U-16B
isoseq3 refine test.5p--3p.tagged.bam 10x_Chromium_3p_primers.fasta test.tagged.refine.bam
isoseq3 correct test.tagged.refine.bam --barcodes 3M-february-2018-REVERSE-COMPLEMENTED.txt test.tagged.refine.corrected.bam
isoseq3 dedup test.tagged.refine.corrected.bam test.tagged.refine.corrected.dedup.bam (optional)
```

We next align the long-reads to the reference genome with pbmm2, as recommended by Iso-Seq.
``` shell
pbmm2 align --preset ISOSEQ --sort test.tagged.refine.corrected.bam ref.genome.fa test.mapped.bam
samtools view -O BAM -F 2052 -h test.mapped.sam |  samtools sort -O BAM -@ 7 -o test.unique.bam -
samtools view -h test.unique.bam | awk '$10 != "*"' |samtools view -bS - > test.flnc.filter.bam
```
