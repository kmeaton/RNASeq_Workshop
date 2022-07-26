#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=03:00:00
#SBATCH --job-name=hisat2

# Load the module
module load gcc/9.3.0
module load hisat2/2.2.1
module load samtools/Jun2022

cd ~

# This first step involves locating splice sites and exons within your genome. It requires that your genome has an annotation, which must be in GTF format. 
# Most genomes that have come out recently will have their annotation information in GFF format, which is slightly different. You can convert between GFF and GTF format using the program gffread. There's great tutorials for this on the internet!

# Identify your splice sites by running:
# hisat2_extract_splice_sites.py A_poly_annotation.gtf > splice_sites.hisat

# Identify your exons by running:
# hisat2_extract_exons.py A_poly_annotation.gtf > exons.hisat

# Now that we have our splice sites and exons, we can use that information to help build a splice site-informed index. This means that the program will take into account splice sites when it is mapping transcripts to the (unspliced) whole genome. 
# Run the following code to build your index:
# The -p flag specifies we want to run this code on 8 threads
# The -ss flag is used to specify the file we just created that contains the known splice sites in our genome
# The --exon flag is used to specify the file we just created that contains the known exons in our genome
# Then, we just list the file containing our genome (in fasta format), and the desired "base name" for our index. The program is going to create several files, but they'll all start with this name. 
# hisat2-build -p 8 -ss splice_sites.hisat --exon exons.hisat Apoly_genome.fa Apoly_genome

# The commented out lines above allow us to build an index, but this can be time consuming, so we've run this code for you ahead of time. 

# Now that we have our index built, let's run HISAT2 to map the reads from each sample to our genome. 
# We'll loop through all of our samples like we did with the trimmomatic script. 
samplelist="Apoly_Dec_1 Apoly_Dec_2 Apoly_Dec_3 Apoly_Dec_4 Apoly_Dec_5 Apoly_Feb_1 Apoly_Feb_2 Apoly_Feb_3 Apoly_Feb_4 Apoly_Feb_5"

# For each sample in our list:
for sample in ${samplelist}
do
	# We'll run hisat2 to align our cleaned reads to the A. polyacanthus genome
	# -p tells the program to run on 8 threads
	# -k tells the program to only search for a maximum of 3 alignments per read. We're only going to continue our analysis with reads that mapped uniquely (i.e., 1 time), so this makes sure the program doesn't spend a bunch of extra time aligning multi-mapping reads
	# -x is the flag that we use to specify the "base name" of our index files
	# -1 is the flag that we use to specify the fastq file containing our mate1 reads
	# -2 is the flag that we use to specify the fastq file containing our mate2 reads
	# After we've created our alignments, we're going to pipe the output of that command to the program samtools, which will convert our output from a .SAM format to a binary .BAM format. The file will encode the same information, but it takes up less space. 
	hisat2 -p 8 -k 3 -x Apoly_genome -1 ${sample}_mate1_paired.fq -2 ${sample}_mate2_paired.fq | samtools view -@ 8 -Sbh > ${sample}.bam
	# We'll use samtools to sort our output alignment (the binary .BAM file), which is necessary to prepare it for the next program we want to use
	samtools sort -@ 8 -o ${sample}_sorted.bam ${sample}.bam
	# Delete the unsorted bam file so it's not taking up unnecessary space
	rm ${sample}.bam
	# Index our sorted bam file, which is necessary for the next program we want to use
	samtools index ${sample}_sorted.bam
done
