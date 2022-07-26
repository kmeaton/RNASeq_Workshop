#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --job-name=trimmomatic

# Load the module
module load java/1.8.0
module load trimmomatic/0.39

# Move into your home directory 
cd ~

# This creates a variable called "samplelist", which contains a list of the names of samples we are going to analyze.
# These names should match the base names of the raw fastq files, removing the "_mate1.fq" or "_mate2.fq"
samplelist="Apoly_Dec_1 Apoly_Dec_2 Apoly_Feb_1 Apoly_Feb_2"

# This allows us to "loop" through every item in the variable $samplelist
# So for each "sample" in our $samplelist, we will run the following command
for sample in ${samplelist}
do
	# This next line of code runs the program trimmomatic
	# java -jar /opt/share/modules/Trimmomatic-0.39/trimmomatic-0.39.jar is what is used to call trimmomatic
	# We specify "PE" after calling trimmomatic because we have paired-end Illumina data (2 files per sample that we sequenced)
	# -threads is used to specify that we want to run this program on 16 cpus
	# -phred33 tells the program that our fastq files use phred 33 encoding. Fastq files contain both sequence and quality information. Phred33 is one method of encoding the quality information for the sequences.
	# After the -phred33 flag, we're specifying the names of the input files and the desired names of our output files. Instead of typing out the entire sample name, we just use ${sample}, which tells the program to use the $sample from our $samplelist
	# For each sample we'll end up with four output files, two paired and two unpaired. 
	# Finally, after specifying the input and output file names, we'll tell the program the trimming parameters we want to use. 
	# ILLUMINACLIP uses known Illumina adapter sequences in a file (TruSeq3-PE-2.fa) and removes them from the reads when detected. 
	# The additional parameters following ILLUMINACLIP are 2:30:10. These can be changed depending on your dataset as necessary, and how stringent you want the filtering to be. 
	# These parameters are used to specify the number of mismatches allowed (2) and the alignment quality thresholds (30 and 10) to declare certain kinds of matches between known adapters and the sequence data we have. 
	# LEADING:4 removes bases with a quality score below 4 from the beginning of the read
	# TRAILING:3 removes bases with a quality score below 3 from the end of the read
	# SLIDINGWINDOW:4:15 trims a read once the average quality within a window of 4 bases falls below an average quality score of 15
	# MINLEN:40 removes all reads with a length below 40 bp
	java -jar /opt/share/modules/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 ${sample}_mate1.fq ${sample}_mate2.fq ${sample}_mate1_paired.fq ${sample}_mate1_unpaired.fq ${sample}_mate2_paired.fq ${sample}_mate2_unpaired.fq ILLUMINACLIP:/opt/share/modules/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:4 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
done

