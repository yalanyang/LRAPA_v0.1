# LRAPA

The pipeline for APA analysis using bulk and single-cell long-read RNA sequencing data.

We would like to construct a data set representing the most comprehensive PAS collection for human brain to date.

Edited in 10/10/2024

**Workflow of LRAPA:**

![](images/APA_workflow.png)

## Installation

You can install LRAPA using command line (linux) by cloning git repository on github:

```         
# Clone the repository 
git clone https://github.com/yalanyang/LRAPA_v0.1.git
cd LRAPA_v0.1
# Create and activate a conda environment
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r

# Install necessary dependencies
conda create -n lrapa -f environment.yml
conda activate lrapa

#in your R session, install these Bioconductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rtracklayer", "Rsamtools", "GenomicAlignments", "Biostrings", "DRIMSeq", "plyranges"))

chmod +x lrapa # Make the script executable
```

## **0. Pre-processing of the long-read CCS reads.**

### **0.1. bulk long-read RNA-seq**

Here, we used one of the three cerebellum samples downloaded from [Pacbio](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/) as test data.

Removal of primers is performed using [lima](https://isoseq.how).

``` shell
lima test.ccs.bam IsoSeq_v2_primers_12.fasta test.fl.bam --isoseq --peek-guess
isoseq3 refine test.fl.IsoSeqX_bc04_5p--IsoSeqX_3p.bam IsoSeq_v2_primers_12.fasta test.refine.bam
```

**Note**: The IsoSeq_v2_primers_12 was download from [Pacbio](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/)

Then we convert the full-length bam file into sam format and run the [flip_reads.py](https://github.com/mortazavilab/ENCODE-references/tree/master) script to orient the reads to the correct strand (since lima\--ccs does not do this, recommend by [Encode long-read pipeline](https://www.encodeproject.org/documents/6d583a1d-d692-4511-b13b-c051822d861c/@@download/attachment/ENCODE%20Long%20Read%20RNA-Seq%20Analysis%20Pipeline%20v3.2%20%28Human%29.pdf)).

``` shell
samtools view -h test.refine.bam > test.refine.sam
python flip_reads.py --f test.refine.sam --o test.fl_flipped.sam
samtools view -bS test.fl_flipped.sam > test.fl_flipped.bam
isoseq3 refine test.fl_flipped.bam PB_adapters.fasta test.flnc.bam
```

**Note**: we didn't trim polyA tails here, because we need to assess whether each read contain a genuine polyA signal based on the length and the base composition of its polyA tail.

We next align the long-reads to the [GRCh38 human reference genome](https://www.gencodegenes.org/human/release_21.html) with Minimap2 using the following parameters:

``` shell
bamtools convert -format fastq -in test.flnc.bam -out test.flnc.fastq
minimap2  -ax splice -uf -C5 $reference/GRCh38.primary_assembly.genome.fa test.flnc.fastq > test.flnc.mapping.sam
samtools view -O BAM -F 2052 -h test.flnc.mapping.sam |  samtools sort -O BAM -@ 7 -o test.flnc.unique.bam -
samtools view -h test.flnc.unique.bam | awk '$10 != "*"' |samtools view -bS - > test.flnc.filter.bam
```

### 0.2. Single-cell long-read RNA-seq

``` shell
lima --isoseq --dump-clips --peek-guess -j 24 test.hifi_reads.bam 10x_Chromium_3p_primers.fasta test.fl.bam
isoseq3 test.fl.5p--3p.bam test.5p--3p.tagged.bam --design T-12U-16B
isoseq3 refine test.5p--3p.tagged.bam 10x_Chromium_3p_primers.fasta test.tagged.refine.bam
isoseq3 correct test.tagged.refine.bam --barcodes 3M-february-2018-REVERSE-COMPLEMENTED.txt test.tagged.refine.corrected.bam
isoseq3 dedup test.tagged.refine.corrected.bam test.tagged.refine.corrected.dedup.bam (optional)
pbmm2 align --preset ISOSEQ --sort test.tagged.refine.corrected.bam ref.genome.fa test.mapped.bam
samtools view -O BAM -F 2052 -h test.mapped.sam |  samtools sort -O BAM -@ 7 -o test.unique.bam -
samtools view -h test.unique.bam | awk '$10 != "*"' |samtools view -bS - > test.flnc.filter.bam
```

## 1. **Filter long-read reads with polyA signals**

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

## **2. PAS identification and annotation**

``` r
lrapa anno -i test.HQ.bam -r $reference/GRCh38.primary_assembly.genome.fa -g $reference/hg38.refGene.gtf -u $reference/hg38.refGene.3utr_merge.bed -o test.PAS.bed
```

**Note:** The hg38.refGene.3utr_merge.bed was generated using the get3UTR_0.1.R script. if you find that this script needs a lot of memory or very slow you may want to split the input bam file by chromosome and run these separately. We do intend to improve this.

**Parameters**

-   `-i, --input`: Input BAM file (required).
-   `-r, --reference`: Reference genome FASTA file (required).
-   `-g, --gtf`: Reference annotation gtf file (required).
-   `-f, --flank_size`: Flank size to search for PAS signal (default: 40).
-   `-u, --utr`: Reference annotation 3'-UTR file (required).
-   `-o, -–output`: Output PAS bed file with header (required).

## **3. PAS count**

### 3.1 Bulk long-read RNA-seq

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

### 3.1 Single-cell long-read RNA-seq

``` shell
lrapa scCount -b scLR.bam -p test.PAS.bed -o RUN3 -w barcode.rev.txt
```

**Parameters**

-   `-i, --input`: Input scRNA-seq BAM file (required).

-   `-p, --pas`: PAS reference TXT file with BED columns followed by annotations (required).

-   `-o, --output`: Output folder for the count matrix (required).

-   `-w, --barcode`: Hitlist barcodes.

## **4. Differential analysis**

### 4.1. With replicates

``` shell
lrapa diff -i test.PAS.count.txt -s sample.txt -c case -n control
```

This module performs differential APA usage analyses between exactly two conditions with two or more replicates based the R package [DRIMSeq](http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html). This is done by testing if the ratio of APA changes between conditions.

If a gene's absolute mean difference of DPAU/gDPAU is \>0.1 and adjusted *P* value is \< 0.01 between two groups, the gene's APA change between two groups will be deemed as significant.

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

### 4.2. No replicates

This module performs differential APA usage analyses between exactly two conditions with no replicate based chisq test. For each test, a n × 2 matrix per PAS site was generated, with the n polyA forming rows and read count in the two conditions forming columns.

We used a Benjamini--Hochberge correction for multiple testing.

If a gene's absolute mean difference of DPAU/gDPAU is \>0.3 and adjusted *P* value is \< 0.01 between two groups, the gene's APA change between two groups will be deemed as significant.

``` shell
lrapa norepdiff -c test.PAS.count.txt -n N1 -g C1
```

**Parameters**

-   `-i, --input`: count file of PAS sites (required).

-   `-c, --group1`: sample 1 for differential analysis (required).

-   `-n, --group2`: sample 2 for differential analysis (required).

    **Note:** The name of sample must same as in the count file.

## 5. Couplings

### 5.1. TSS-PAS coupling

![](images/TSS-polyA.png)

We counted 5ʹ-3ʹ isoforms using GenomicFeatures. Each Pacbio cDNA read was assigned to a TSS in a window of 50 nt and to a PAS. Only the reads that mapped to both features were retained and considered full-length reads. Counts were summarized in 5ʹ-3ʹ isoforms, resulting in counts for each 5ʹ-3ʹ combination.

``` shell
lrapa couplingTSS -i test.mapping.bam -g hg38.refGene.gtf -p test.PAS.bed -o TSS-PAS.coordination.txt
```

TSS: Transcriptional start site

**Parameters**

-   `-i, --input`: Input full-length BAM file (required).
-   `-g, --gtf`: Reference annotation gtf file (required).
-   `-p, --pas`: Reference pas bed (required).
-   `-o, --output`: Output TSS-exon_coordination file.

Reads were filtered to retain only full-length reads: The script will generated a file (full-length.reads.txt) contained the qname of the full-length reads (TSS+TES). Then we get extract the full-length reads from the bam file using samtools, which will be used for PAS-exon coupling analysis.

``` shell
samtools view -H test.mapping.bam > test.header
samtools view -N full-length.reads.txt -o test.full_length.sam test.mapping.bam
cat test.header test.full_length.sam  > test.full_length.header.sam 
samtools view -bS test.full_length.header.sam > test.full_length.bam
```

### 2. Exon-PAS coupling

![](images/Exon-polyA.png)

We quantify the regulatory links between exons and 3ʹ ends. Given that every read represents a full-length transcript, we assessed all features of each read to quantify the frequency of co-occurrence between features using χ2 test.

At this step, we first extract skipped exons from the reference gtf annotation file or assemble transcriptome gtf provided by the users. The exons that are overlapped with other exons and 3'-UTR are removed.

Using this read to feature assignment, we counted the number of reads assigned to a given polyA site and divided them into reads including a particular exon or skipping the exon. Testing for exon--end site coordination were performed using a χ2 test. For each test, a n × 2 matrix per SE was generated, with the n polyA forming rows and inclusion and exclusion counts forming columns. Finally, we used a Benjamini--Hochberge correction for multiple testing and reported the FDR value.

``` shell
lrapa couplingExon -i test.full_length.bam -g hg38.refGene.gtf -p test.PAS.bed -o PAS-exon.coordination.txt
```

**Options**

-   `-i, --input`: Input full-length BAM file, generated by .

-   `-g, --gtf`: Reference annotation gtf file (required).

-   `-p, --pas`: Reference pas bed.

-   `-o, --output`: Output TSS-exon_coordination file.

## File conversion scripts

get3UTR_0.1.R : get the 3'UTR reference for PAS analysis

## **Additional programs**

**References:**

[Alfonso-Gonzalez C, Legnini I, Holec S, et al. Sites of transcription initiation drive mRNA isoform selection[J]. Cell, 2023, 186(11): 2438-2455. e22.](https://www.cell.com/cell/fulltext/S0092-8674(23)00408-7)

[Hardwick S A, Hu W, Joglekar A, et al. Single-nuclei isoform RNA sequencing unlocks barcoded exon connectivity in frozen brain tissue[J]. Nature biotechnology, 2022, 40(7): 1082-1092.](https://www.nature.com/articles/s41587-022-01231-3)
