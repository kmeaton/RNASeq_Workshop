#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=03:00:00
#SBATCH --job-name=stringtie

module load stringtie/2.2.1

# Move into your home directory
cd ~/

# Make a subdirectory where we are going to output the results file
mkdir stringtie_results

# Now we will run the program stringtie, which will use our alignments (the sorted .bam files) to estimate expression levels for each gene in our genome

# The -p flag tells the program to run on 8 threads
# -e tells stringtie to operate in expression estimation mode, limiting the processing of read alignments to estimating the coverage of transcripts given in the -G option
# -G allows you to specify a gff file containing annotation information for the genome you mapped your reads to
# -B specifies the output format of tables containing coverage data for the reference transcripts in the GFF file

samplelist="Apoly_Dec_1 Apoly_Dec_2 Apoly_Feb_1 Apoly_Feb_2"

for sample in ${samplelist}
do
	mkdir stringtie_results/${sample}
	stringtie -p 8 -e -B -G genomic.gff -o ~/stringtie_results/${sample}/${sample}.gtf ${sample}_sorted.bam
done
