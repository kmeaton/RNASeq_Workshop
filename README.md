# RNASeq Workshop
This repository contains a walkthrough of how to analyze RNA-Seq data, using a sample dataset from Bernal et al. 2020 (*Science Advances*). This repository will be used in the RNA-Seq data analysis workshop held at UT Chattanooga from August 10-12, 2022.

## Contents
1. [Getting set up](https://github.com/kmeaton/RNASeq_Workshop#getting-set-up)
2. [Cleaning and mapping reads](https://github.com/kmeaton/RNASeq_Workshop#day-1-cleaning-and-mapping-reads)
3. [Analyzing expression patterns](https://github.com/kmeaton/RNASeq_Workshop#day-2-testing-for-differential-gene-expression)

## Getting set up

If you don't already have an account on [Github](https://github.com), make one now.

1. Log in to your [Github](https://github.com) account.

2. Set up a Personal Access Token so that you can access your Github account from the command line. Go to the top right-hand corner of the page which has a circle (containing your profile picture). Click on this icon, then click "Settings". From your settings page, scroll down the menu of options on the left hand side of the screen. Click on "Developer settings". Then, again in the menu on the left side of the page, click on "Personal access tokens". Click "Generate new token". Under "Note", name this personal access token "RNASeq-Workshop". Set the key to expire in 30 days (the default. Under "Select scopes", click the first box, next to "repo". Then scroll all the way down the page and click the green "Generate token" button. This will essentially create a temporary password for you to access your Github account from the command line. Copy the personal access token you have just generated (it should just be a string of random letters and numbers). Make sure you save it somewhere, because you won't be able to go back and re-copy it to your clipboard later! I recommend copy/pasting it into the notes app on your computer, or into a Word document you have saved. Don't close out of the page until you've got the token saved. 

3. Once your token is saved, navigate back to this page. Fork this repository by clicking the "Fork" button on the upper right of this page. After a few seconds, you should be looking at your own copy of this repository in your own Github account. 

4. Click the green "Code" button at the upper right of this page. Click the tab that says HTTPS, then copy the link that's shown below it. 

5. Log in to your UTC computing cluster account by typing the following code into the terminal, substituting your UTC username in where it says [user]. You'll be prompted to enter a password, which you'll type right into the terminal.
```shell
ssh [user]@epyc.simcenter.utc.edu
```

6. Once you're logged in, in your home directory, type the following to clone into the repository. Make sure you're cloning into __your__ fork of the repository, not my original one.
```shell
git clone the-url-you-copied-in-step-3
```

7. It will prompt you to enter your Github username. Type it right into the terminal. Then, it will prompt you to enter your password. Copy the Personal Access Token you generated in Step 2, and paste it into the terminal, then hit enter. This will create a copy of this repository in your account on the Epyc cluster. 

8. Next, move into this directory:
```shell
cd RNASeq_Workshop
```

9. At this point, you should be in your own local copy of the repository, which contains all the scripts you'll need to edit and run to analyze our practice dataset. 

10. Finally, let's make the scripts that you're going to run executable. From inside the ```RNASeq_Workshop``` folder, type the following:

```shell
chmod +x *.sh
```

**Congrats!** You're all set up to process some RNASeq data!

## Day 1: Cleaning and mapping reads

Log in to your account on the UTC Computing Cluster, substituting your UTC username in where it says [user]. 

```shell
ssh [user]@epyc.simcenter.utc.edu
```

You'll be prompted to enter your password, and then you'll be logged in to the cluster. Move into the directory you cloned yesterday, which contains all the scripts you'll need to work with today.

```shell
cd RNASeq_Workshop
```

### Assessing raw read quality

For this tutorial, we are going to use a sample RNA-Seq dataset from [this paper](https://www.science.org/doi/full/10.1126/sciadv.aay3423). The authors examined how gene expression changed in several species of fish as they were exposed to unseasonably warm temperatures over several months. We will use RNA-Seq data from one of the species (the spiny chromis damselfish, *Acanthochromis polyacanthus*) and compare gene expression between two months (December, when temperatures were relatively normal, and February, when temperatures were far above average). 

In the shared class directory, there are 8 files with the extension ```.fq```. These fastq files contain "raw" RNA-Seq reads for 4 samples - 2 from December and 2 from February. Each sample will have data in two different files: one will have the extension ```_mate1.fq```, and the other will be ```_mate2.fq```. This is because these samples were sequenced using paired-end Illumina sequencing, which generates two sequences per input molecule of RNA (one sequence in the "forward" direction, and one in the "reverse" direction). These paired sequences are stored in two different files.

You'll need to copy these files from the class directory to your home directory before you can start working on them. Type the following code in your terminal:

```shell
cp /scr/south_east_comp/*.fq ~/
cd ~
ls
```

You should now see that in your home directory, you have your own copies of the raw data files. 

Before we can do any analysis on our raw data, we have to clean them to remove any low-quality reads or contamination. We'll start by examining the quality of the sequences, using the program [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 

First, create a new directory for our FastQC results to go into:
```shell
mkdir ~/raw_reports
```

Then, move into the RNASeq_Workshop folder that you cloned from Github. This folder contains all the scripts you'll need to run today. Examine the code in the script ```fastqc.sh```, and then run the script by typing the following:
```shell
# First examine the code in the script
cat fastqc.sh
# Once you understand what the code is doing, run the script
sbatch fastqc.sh
```

This should take just a couple of minutes to run. When the script is done running, take a look at the output. You'll have to download the reports to your local machine. Open a terminal on your __local machine__ and type the following:
```shell
scp -r [user]@epyc.simcenter.utc.edu:~/raw_reports/*.html ~/Desktop/
```

You'll be prompted to enter your UTC password here, and then the reports will be copied from your fastqc run to the Desktop on your local computer. Click on the files that you just copied, open them up, and see what the sequencing data looks like. Is it high quality? How can you tell? 

### Removing adapters and low-quality sequences

Now that we know what our raw reads look like, we should trim the sequencing adapters and remove any low-quality reads. We can do this using the program Trimmomatic. Move into the RNASeq_Workshop folder and examine the Trimmomatic script we are going to run. For more information on the quality trimming parameters, examine [the manual for trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf). Once you understand what the code is doing, run the script by typing the following:

```shell
# Examine the code in the script
cat trimmomatic.sh
# Run the script
sbatch trimmomatic.sh
```

This should only take a few minutes. This program takes paired-end sequencing data for each of our samples, trims any sequencing adapters that are present in the reads, and then removes any low-quality sequences from the dataset. You'll end up with four files for each sample - here's an example of what it'll look like for the sample ```Apoly_Dec_1```:

```
Apoly_Dec_1_mate1_paired.fq
Apoly_Dec_1_mate1_unpaired.fq
Apoly_Dec_1_mate2_paired.fq
Apoly_Dec_1_mate2_unpaired.fq
```

We will continue our analyses only using the "paired" files. These contain sequences where both members of a "pair" of reads have passed the quality filtering step. Check out the sizes of the fastq files containing the paired and unpaired reads by running the following commands:

```shell
ls -lh *_paired.fq
ls -lh *_unpaired.fq
```

The fourth column in the output from this command will show the approximate size of the file. Since the majority of our reads were of high quality, the "paired" files (where both reads in a mate pair passed the quality filtering step) should be several times larger than the "unpaired" files. 

### Checking the quality of our trimmed and filtered reads

Let's make sure that the filtering steps worked well, and that the quality of our sequences increased after trimming and filtering. To do this, we'll run FastQC again, but this time on the ```_paired.fq``` files, which contain the sequences we will proceed with for our analyses.

We'll make a new directory for these reports.

```shell
mkdir filtered_reports
```

Then, modify the ```fastqc.sh``` script in your folder so that it will run on your filtered files and send the output to the filtered_reports folder we just created, and save it as a new file called ```fastqc_filtered.sh```. You can do this in your favorite text editor, like nano or vim. 

Run your new script, ```fastqc_filtered.sh```, and then download the reports to your local machine like we did before. Check out the sequence quality now - how has it improved?

### Mapping our trimmed reads to a reference

Now that we have high-quality, filtered reads for each of our samples, we need to map them to a reference genome. This allows us to quantify the number of sequences in our dataset that originated from each gene in our organism's genome. For our analyses, we will use the published genome of the spiny chromis damselfish, *A. polyacanthus* (publicly available on NCBI at accession number: GCF_002109545.1). 

We will use a program called [HISAT2](http://daehwankimlab.github.io/hisat2/manual/) to map the RNA-Seq reads to the *A. polyacanthus* genome. HISAT2 requires that you first build an index of the reference genome, which it then uses during the mapping process. Indexing the reference genome takes about an hour, so we have provided the indexed genome for you. The code that we used to index the genome is included in the script ```hisat2.sh```, so that you can see how it was done, but it has been commented out so that you don't have to run it. 

Instead of building our index, we'll just copy it from the shared class folder. 

```shell
cp /scr/south_east_comp/*.ht2 ~/
```

Once your index is copied from the shared folder, examine the ```hisat2.sh``` script in your RNASeq_Workshop folder. Once you understand what it is doing, run it.

```shell
# Examine the script
cat hisat2.sh
# Run the script
sbatch hisat2.sh
```

Formatting the data into sorted, indexed ```.bam``` files should take 10-15 minutes. 

When this script is done running, take a look at the SLURM output file (it will be called ```slurm-[XXXXX].out```, where [XXXXX] is the job number). Examine the overall mapping percentage for each of your four samples. Does it look like we have a good reference? How can you tell?

### Generating read count data

For each sample, we have mapped our RNA-Seq reads to our reference genome. We now need to generate a matrix of counts that correspond to the expression levels of each gene. We can do this using the program [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual). 

StringTie will look at our alignment and use the genome annotation information that we provide in ```genomic.gff``` to count the number of RNA-Seq reads that have mapped to each gene and transcript in our *A. polyacanthus* genome. You'll need to copy this ```genomic.gff``` file from the shared directory to your home directory by doing the following:

```shell
cp /scr/south_east_comp/genomic.gff ~/
```

Then, once you've copied the annotation file to your home directory, examine the ```stringtie.sh``` script in your RNASeq_Workshop folder. Once you understand what it is doing, run it.

```shell
# Examine the script
cat stringtie.sh
# Run the script
sbatch stringtie.sh
```

### Formatting read count data for DESeq2

Finally, you'll need to format your read count data to be appropriately read in to DESeq2, the R package we will use for our differential gene expression tests. The developers of StringTie provide a Python script called ```prepDE.py``` on their website for this purpose. I've copied the script to our RNASeq_Workshop folder, but you can also download it directly from the source [here](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py). 

```shell
# Start by moving into your RNASeq_Workshop folder
cd ~/RNASeq_Workshop
# Run the prepDE.py script from the command line
# The -l flag allows us to specify that the average read length in our alignment was 151 bp
# The -i flag tells the script where to look for the input files. Our stringtie outputs for each sample were written to a sample-specific folder in ~/stringtie_results. 
python prepDE.py -l 151 -i ~/stringtie_results/
```

You should now have two output files from this command: one called ```gene_count_matrix.csv``` and one called ```transcript_count_matrix.csv```. These files contain the number of reads that mapped to each gene (or transcript) in the *A. polyacanthus* genome, for each sample. The files are formatted like this:

| Gene | Apoly_Dec_1 | Apoly_Dec_2 | Apoly_Feb_1 | Apoly_Feb_2 |
| --- | --- | --- | --- | --- |
| gene-1 | 12 | 15 | 40 | 38 | 
| gene-2 | 705 | 813 | 30 | 55 | 
| gene-3 | 0 | 2 | 44 | 17 |

You can see that each sample has its own column, and each gene has its own row. For our differential gene expression analyses tomorrow, we will use the ```gene_count_matrix.csv``` file, because we are interested in differential gene expression, and not particularly interested in differential expression patterns of transcripts. If you were doing a study where transcript-specific differences in expression patterns were important, you might want to use the ```transcript_count_matrix.csv``` file instead. 

## Day 2: Testing for differential gene expression

Today, we are going to use a gene count matrix to test for statistically significant differentially expressed genes between fish in two different treatment groups. We prepared one of these gene count matrices yesterday, for four samples. Because it's harder to detect significant differences between treatments when you have fewer samples per treatment, today we'll use a gene count matrix that has ten samples (five from each treatment). This matrix is on the class GitHub page, it's called ```gene_count_matrix_full.csv```. Download this file to your **local** machine. We won't need to do anything on the computing cluster today!

Download this file by scrolling to the top of the page, clicking on the ```gene_count_matrix_full.csv``` file, and then clicking the "Download" button towards the right hand side of the screen. Do the same thing with the script ```DESeq2.R```, as well as the table that contains information about each of our samples, ```sample_info.csv```. These are the three files we'll be working with today. 

Open up RStudio on your computer, and open the ```DESeq2.R``` script you just downloaded. Read through the script and follow along with the instructions in the comments to analyze the expression data.

## Acknowledgements
This workshop was made possible by funding provided to Fernando Alda from the University of Tennessee at Chattanooga. 

The sample datasets used in this analysis are publicly available on NCBI. The version of the *A. polyacanthus* genome used is available at accession number GCF_002109545.1. The raw RNA-Seq reads are available under NCBI BioProject Number PRJNA489934 and SRA accession number SRP160415. We thank the following authors of the study that generated this dataset: Moises Bernal, Celia Schunter, Robert Lehmann, Damien Lightfoot, Bridie Allan, Heather Veilleux, Jodie Rummer, Philip Munday, and Timothy Ravasi, as they have graciously allowed us to use their data (Bernal et al. 2020, *Science Advances* 6(12): eaay3423). 

The formatting of this tutorial, as well as the "Getting set up" portion of this page, borrow heavily from the Scripting for Biologists python tutorials written by [Jamie Oaks](http://phyletica.org) available [here](https://github.com/joaks1). 

Portions of the analysis pipeline followed here (especially the code for StringTie, the code used to prep read counts for DESeq2, and portions of the R code and analyses in DESeq2) were based on code and activities from [Tonia Schwartz's](https://www.schwartzlab-ecoevolutionarygenomics.org/) Functional Genomics class at Auburn University.
