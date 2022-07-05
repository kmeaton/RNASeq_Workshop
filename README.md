# RNASeq Workshop
This repository contains a walkthrough of how to analyze RNA-Seq data, using a sample dataset from Bernal et al. 2020 (*Science Advances*). This repository will be used in the RNA-Seq data analysis workshop held at UT Chattanooga from August 10-12, 2022.

## Contents
1. [Getting set up](https://github.com/kmeaton/RNASeq_Workshop#getting-set-up)
2. [Cleaning and mapping reads](https://github.com/kmeaton/RNASeq_Workshop#day-1-cleaning-and-mapping-reads)
3. [Counting reads and analyzing expression patterns](https://github.com/kmeaton/RNASeq_Workshop#day-2-generating-read-counts-and-testing-for-differential-gene-expression)

## Getting set up


## Day 1: Cleaning and mapping reads

Log in to your account on the UTC Computing Cluster. 

```shell
ssh [user]@[UTC address]
```

You'll be prompted to enter your password, and then you'll be logged in to the cluster. 

In your home directory, you should have four fastq files with the extension ".fq". These are our raw sequences. Before we can do any analysis on them, we have to clean them to remove any low-quality reads or contamination. We'll start by examining the quality of the raw sequences, using the program FastQC. 

First, create a new directory for our FastQC results to go into:
```shell
mkdir raw_reports
```

Then, run the script ```fastqc_raw.sh``` in the directory with your fastq files. This script should take around 20 minutes to run on your 4 files. 

When the script is done running, take a look at the output. You'll have to download the reports to your local machine. Open a terminal on your __local machine__ and type the following, substituting in your UTC username where it says [user]:
```shell
scp -r [user]@[UTC address]:/home/[user]/raw_reports/*.html ~/Desktop/
```

This should copy the output files from your FastQC run to your Desktop on your local computer. Click on the files that you just copied, open them up, and see what the sequencing data looks like. Is it high quality? How can you tell? 



## Day 2: Generating read counts and testing for differential gene expression


## Acknowledgements
This workshop was made possible by funding provided to Fernando Alda from the University of Tennessee at Chattanooga. 

The sample datasets used in this analysis are publicly available on NCBI. The version of the *A. polyacanthus* genome used is available at accession number GCF_002109545.1. The raw RNA-Seq reads are available under NCBI BioProject Number PRJNA489934 and SRA accession number SRP160415. We thank the following authors of the study that generated this dataset: Moises Bernal, Celia Schunter, Robert Lehmann, Damien Lightfoot, Bridie Allan, Heather Veilleux, Jodie Rummer, Philip Munday, and Timothy Ravasi, as they have graciously allowed us to use their data (Bernal et al. 2020, *Science Advances* 6(12): eaay3423). 

