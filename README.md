# RNASeq Workshop
This repository contains a walkthrough of how to analyze RNA-Seq data, using a sample dataset from Bernal et al. 2020 (*Science Advances*). This repository will be used in the RNA-Seq data analysis workshop held at UT Chattanooga from August 10-12, 2022.

## Contents
1. [Getting set up](https://github.com/kmeaton/RNASeq_Workshop#getting-set-up)
2. [Cleaning and mapping reads](https://github.com/kmeaton/RNASeq_Workshop#day-1-cleaning-and-mapping-reads)
3. [Counting reads and analyzing expression patterns](https://github.com/kmeaton/RNASeq_Workshop#day-2-generating-read-counts-and-testing-for-differential-gene-expression)

## Getting set up

If you don't already have an account on [Github](https://github.com), make one now.

1. Log in to your [Github](https://github.com) account.

2. Fork this repository by clicking the "Fork" button on the upper right of this page. 
..After a few seconds, you should be looking at your own copy of this repository in your own Github account. 

3. Click the green "Code" button at the upper right of this page. Click the tab that says HTTPS, then copy the link that's shown below it. 

4. Log in to your UTC computing cluster account by typing the following code into the terminal, substituting your UTC username in where it says [user]. You'll be prompted to enter a password, which you'll type right into the terminal.
```shell
ssh [user]@[UTC address]
```

5. Once you're logged in, in your home directory, type the following to clone into the repository. Make sure you're cloning into __your__ fork of the repository, not my original one.
```shell
git clone the-url-you-copied-in-step-3
```

6. Next, move into this directory:
```shell
cd RNASeq_Workshop
```

7. At this point, you should be in your own local copy of the repository, which contains all the scripts you'll need to edit and run to analyze our practice dataset. 

## Day 1: Cleaning and mapping reads

Log in to your account on the UTC Computing Cluster, substituting your UTC username in where it says [user]. 

```shell
ssh [user]@[UTC address]
```

You'll be prompted to enter your password, and then you'll be logged in to the cluster. Move into the directory you cloned yesterday, which contains all the scripts you'll need to work with today.

```shell
cd RNASeq_Workshop
```

### Assessing raw read quality

In your home directory (MAYBE? MIGHT HAVE TO COPY THESE IN? CAN WE HAVE A CLASS DIRECTORY?), you should have four fastq files with the extension ".fq". These are our raw sequences. Before we can do any analysis on them, we have to clean them to remove any low-quality reads or contamination. We'll start by examining the quality of the raw sequences, using the program FastQC. 

First, create a new directory for our FastQC results to go into:
```shell
mkdir raw_reports
```

Then, run the script ```fastqc_raw.sh``` in the directory with your fastq files. This should take about 15-20 minutes. 
```shell
bash fastqc_raw.sh
```

When the script is done running, take a look at the output. You'll have to download the reports to your local machine. Open a terminal on your __local machine__ and type the following, substituting in your UTC username where it says [user]:
```shell
scp -r [user]@[UTC address]:/home/[user]/raw_reports/*.html ~/Desktop/
```

This should copy the output files from your FastQC run to your Desktop on your local computer. Click on the files that you just copied, open them up, and see what the sequencing data looks like. Is it high quality? How can you tell? 

### Removing adapters and low-quality sequences

Now that we know what our raw reads look like, we should trim the sequencing adapters and remove any low-quality reads. We can do this using the program Trimmomatic. Examine the Trimmomatic script we are going to run by typing the following:

```shell
cd RNASeq_Workshop
less trimmomatic.sh
```

## Day 2: Generating read counts and testing for differential gene expression


## Acknowledgements
This workshop was made possible by funding provided to Fernando Alda from the University of Tennessee at Chattanooga. 

The sample datasets used in this analysis are publicly available on NCBI. The version of the *A. polyacanthus* genome used is available at accession number GCF_002109545.1. The raw RNA-Seq reads are available under NCBI BioProject Number PRJNA489934 and SRA accession number SRP160415. We thank the following authors of the study that generated this dataset: Moises Bernal, Celia Schunter, Robert Lehmann, Damien Lightfoot, Bridie Allan, Heather Veilleux, Jodie Rummer, Philip Munday, and Timothy Ravasi, as they have graciously allowed us to use their data (Bernal et al. 2020, *Science Advances* 6(12): eaay3423). 

The formatting of this tutorial, as well as the "Getting set up" portion of this page, borrow heavily from the Scripting for Biologists python tutorials written by [Jamie Oaks](http://phyletica.org) available [here](https://github.com/joaks1). 
