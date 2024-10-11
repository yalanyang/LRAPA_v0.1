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

The complete LRAP manual is available via [github page](https://yalanyang.github.io/LRAPA_v0.1).

More details could be found in the [tutorial](https://yalanyang.github.io/LRAPA_v0.1/).

If you have any questions about LRAPA, please directly contact Yalan Yang (yangyalan\@uchicago.edu).
