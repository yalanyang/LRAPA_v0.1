---
layout: default
title: Installation
nav_order: 2
---

# **Installation**

You can install LRAPA using command line (Linux or macOS) by cloning the repository on github:

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

# Make the script executable
chmod +x lrapa 
```
