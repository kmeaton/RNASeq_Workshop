# DESeq2 script for analysis of differentially expressed genes
# Written by: KM Eaton, Auburn University, 2022
# Big thanks to Dr. Tonia Schwartz and her Functional Genomics course at Auburn University for providing some of the code upon which this tutorial is based!
# Manual for DESeq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Begin by setting your working directory to wherever you downloaded the gene count matrix and the sample info table.
# This should be the same folder that this script is in!
setwd("~/Desktop/")

# Install the DESeq2 package if you haven't already. If it's already installed, skip ahead to loading the library.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Load the DESeq2 library.
# Unlike installing the package, you'll need to load the library every time you run this script, not just the first time.
library(DESeq2)

#################     Input data      #################

# After the packages we need are installed and loaded, we'll read in our data.

# We'll make an object called countdata - this will be our gene count matrix that has read counts for each gene and each sample.
countdata<-as.matrix(read.csv("gene_count_matrix_full.csv", row.names="gene_id"))
# Examine the first few lines of the matrix to make sure that it looks right:
head(countdata)

# Next, we'll make an object called coldata - this will be the file that contains information on our samples and which treatment group each sample belongs to
coldata<-read.csv("sample_info.csv", header = TRUE, row.names = 1)
# Examine the first few lines of the table to make sure that it looks right:
head(coldata)

# Now make sure that all sample IDs in your coldata table are also in your countdata table
all(rownames(coldata) %in% colnames(countdata))
# If this does not return TRUE, go back and make sure that every sample in your coldata table (the one with the sample and treatment information) was included in your count matrix
all(rownames(coldata) == colnames(countdata))
# If this does not return TRUE, the ORDER of the samples in your coldata table is incorrect. Reorder the samples in this coldata table so that they are in the same order (from top to bottom) as the samples in your countdata matrix (from left to right)

# If both of the previous lines of code returned TRUE, proceed to creating a DESeq dataset. 
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design =~ treatment)
# When creating a DESeq dataset, you always must specify your countData (in the form of a matrix where the rows are genes and the columns are samples, and each box contains the number of reads that mapped to a particular gene in a particular sample)
# You must also always specify your colData, in the form of a table that contains metadata information about each sample. At the very minimum, it must contain two columns - one with the sample name, and one with the "condition" or "treatment" applied to that sample
# Third, you must always specify a design formula. If you've only got one independent variable (i.e., one thing that changed between your groups of interest), this should be straightforward. 
# In our case, we will just specify design =~treatment. Every variable on the right side of the equals sign in your design formula must be the same as the name of a column in your coldata table (i.e., an attribute that each of your sequenced samples has)

# Take a look at your DESeq dataset
dds
# The second line of this output should say "dim:" followed by two numbers. 
# These are the "dimensions" of our dataset. What do you think each of these numbers represents? 

#################     Differential expression analysis      #################

# To run the statistical analysis on our dataset, we use the command DESeq()
# This analyzes how gene expression has changed between the control and experimental groups. 
# For this dataset, our control group is going to be the "December" samples, since they were experiencing cool temperatures that are relatively normal.
# Our experimental group will be the "February" samples, because they were experiencing heatwave conditions with temperatures well above the seasonal average. 
# First, we need to tell DESeq that December is our "reference" or "control" group. 
dds$treatment<-relevel(dds$treatment, ref = "December")

# Now that we've specified which group is our reference, we can run DESeq().
dds<-DESeq(dds)

# This should take a couple of seconds to run, so don't worry if it doesn't finish immediately!

# To look at the results, we'll use the function results()
# The results function takes a couple of different arguments. Importantly, you'll want to specify alpha = 0.05. The default value of alpha is 0.1, but you want to use a value of alpha that matches the p-value cutoff you're going to impose to detect significantly differentially expressed genes. 
# Since our p-value threshold will be 0.05, alpha should also be 0.05
# We'll also use the contrast argument to tell the program which two groups to compare in our dataset. 
# Since we only have two different treatment groups, this isn't super important here, but it becomes more important when you have >2 levels of your independent variable.
# This contrast argument takes three values - the variable in the design formula that you want to compare results by (in our case, treatment), the "experimental" group, and the "control" group. 
res<-results(dds, alpha = 0.05, contrast = c("treatment", "February", "December"))
# Take a look at your results - this will just output several lines of a big table.
res
# Let's look at the summary information for this table:
summary(res)
# How many total expressed genes were we able to detect in our dataset?
# LFC > 0 (up) represents the number of genes that were upreglated (i.e., more highly expressed) in the February samples as compared to the December samples
# LFC < 0 (down) represents the number of genes that were downregulated (i.e., less highly expressed) in the February samples as compared to the December samples

# How about we subset our results to only include signficantly differentially expressed genes between the two treatment groups?
res_sigs<-subset(res, padj < 0.05)
# Why do we have to use p-adj rather than the raw p-value?
# Take a new look at the results - this will again just output a few lines of a big table
res_sigs
# This table contains just the genes that were significantly differentially expressed between the December and February groups (p-adj < 0.05)
# What does each of the columns in this table represent?

#################     Generating expression plots      #################

# This next function, plotMA() creates a scatterplot of log2fold changes in expression level (i.e., how strongly up- or down-regulated a particular gene is) vs. the mean of normalized counts (i.e., how highly expressed a gene is on average)
DESeq2::plotMA(res, main = "Log2-Fold Change vs. Base Mean", ylim=c(-5,5))
# Each point on the plot represents a gene. The significantly differentially expressed genes are colored blue. Triangles pointing up or down fall outside of the limits of the plot.
# Can you generate this plot but make it only show the significant genes?

# Some types of plots require that you transform your expression data to normalize the variance
# This transformation removes the experiment-wide trend of variance in count data 
# You can do this using a variance stabilizing transformation (vst), or you can use a regularized logarithm (rlog)
# These transformations essentially do the same thing, but vst is generally a bit faster, so that's what we'll use here
vsd<-varianceStabilizingTransformation(dds)

# We can generate a heatmap of our count data, which will cluster samples together based on their similarity in expression patterns:
# First, we must install the package "pheatmap" (we only need to run this next line once - after pheatmap is installed, you shouldn't ever have to install it again)
install.packages("pheatmap")
# Now load the library
library(pheatmap)

# We'll also install the package "RColorBrewer" which will allow us to color our heatmap nicely!
install.packages("RColorBrewer")
# Load the library
library(RColorBrewer)

# Next, we can calculate the "sample-to-sample distances" (how different are any two samples from one another, in terms of their gene expression patterns)
sampleDists <- dist(t(assay(vsd)))
# Turn this vector (i.e., list of numbers) into a matrix
sampleDistMatrix <- as.matrix(sampleDists)
# Make the row names of your matrix equal to the values in the "treatment" column of your vsd object
rownames(sampleDistMatrix) <- paste(vsd$treatment)
colnames(sampleDistMatrix) <- NULL # remove column names in your matrix to make it less cluttered
# Generate a pretty color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
# Generate your heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists,
         col=colors)

# Examine the heatmap we just generated. The scale on the right hand side of the figure shows the "sample distance". 
# The numerical values on this scale are relatively meaningless, but what's important to note is that larger values mean samples were more distant from one another (i.e., had less similar patterns of gene expression). 
# In the color palette we've chosen, samples that are more similar (i.e., have lower sample-to-sample distances) are darker blue, whereas samples that are more different are lighter blue/white. 
# You'll note that there are two large squares of relatively darker blue on our heatmap. Why is this? Is this something you would expect, given what you know about the dataset?

# Next, we can generate a principal component analysis (PCA) of our gene expression data. 
# This PCA also clusters samples based on similarity in expression patterns. 

# We'll use a function built-in to DESeq2 that's able to handle this format of data:
plotPCA(vsd, intgroup=c("treatment"))
# What does this PCA tell us about our samples? Is this what you would expect to see?

#################     Writing results      #################

# Let's save our gene expression results into two different files: one that contains the expression results for all of the expressed genes, and one that contains the expression results for just the significantly differentially expressed genes
# First, we'll do all the genes
write.csv(res, file = "gene_expression_results_full.csv")
# Now, just the significantly differentially expressed genes
write.csv(res_sigs, file = "gene_expression_results_sigs.csv")

#################     Exploring results      #################

# Open this output CSV file in excel. Sort the data by log2-fold change. Try googling some of the most highly downregulated (genes with the most negative log2-fold change)
# and most highly upregulated (genes with the most positive log2-fold change) genes. Are these genes that you expect would be up- or down-regulated under heat stress? 
# Don't try googling genes that have names like "LOC0000000" - these are unannotated genes with no functional information. We know they are genes in the A. polyacanthus genome, but have no idea what they are or what they do.
# Now try searching for genes that we would definitely expect to be up- or down-regulated under heat stress. (Hint: Try searching "HSP".) Are these genes in our list of significant DEGs?
# Are they significant in the way you'd expect?

#################     WGCNA      #################

# WGCNA (weighted gene co-expression network analysis) is an R package that can be used to identify
# co-expressed modules of genes (i.e., genes that are often expressed together in your dataset). 
# This portion of our script is based on a tutorial available here: https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#1_Purpose_of_this_analysis
# Also based on code from Ally Swank, graduate student at Auburn University

# Use the code below to install and load the required packages for this analysis (there's a lot!)

BiocManager::install("impute")
library(impute)

BiocManager::install("GO.db")
library(GO.db)

BiocManager::install("preprocessCore")
library(preprocessCore)

install.packages("WGCNA")
library(WGCNA)

install.packages("readr")
library(readr)

install.packages("ggforce")
library(ggforce)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

install.packages("dplyr")
library(dplyr)

if (!("magrittr" %in% installed.packages())) {
  install.packages("magrittr")
}
library(magrittr)

if (!("limma" %in% installed.packages())) {
  BiocManager::install("limma")
}
library(limma)

if (!("ggplot2" %in% installed.packages())) {
  install.packages("ggplot2")
}
library(ggplot2)

# Add a column to your "coldata" table (which contains the treatment information for each fish) stating the fish's ID
# This will be the same as the rowname in the coldata table
coldata$fishID <- rownames(coldata)

# Retrieve the normalized expression data from the "vsd" object (which is in "DESeqDataSet" format) and convert it into a matrix (using the assay function).
# This gives you a matrix where the samples are columns, and the genes are rows. Transpose this using the function t(), to switch the rows/columns (i.e., make it so samples are rows and genes are columns)
normalized_counts <- assay(vsd) %>%
  t()

# To identify which genes are co-expressed similarly across samples, WGCNA first creates a weighted network to define which genes are most similar in expression.
# The measure of similarity it uses is based on a correlation matrix, but it requires the definition of a threshold value, which in turn depends on a "power" parameter that defines the exponent used when transforming the correlation values. 
# The choice of power parameter will affect the number of modules identified. 
# We can use this pickSoftThreshold function to identify a good choice for our power parameter. 
# This code might take a little while to run (10 minutes or so)
sft <- pickSoftThreshold(normalized_counts, dataIsExpr = T, networkType = "signed")

# Calculate the signed R^2, a measure of the model fit. We will plot this value to figure out what our "power" threshold should be. 
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

# Now we'll plot the "model fit" (i.e., the signed R^2) vs. the power threshold to decide which power threshold gives us the best model fit. 
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) + # This line passes our data to the program ggplot
  geom_point() + geom_text(nudge_y = 0.1) + # This tells it that we want it to be a dotplot, and that we want each point to have a label right next to it
  geom_hline(yintercept = 0.80, col = "red") + # This plots a horizontal red line on the same axes, according to the formula y = 0.8 (which is the signed R^2 threshold that the program recommends)
  ylim(c(min(sft_df$model_fit), 1.05)) + # This sets the bounds of our y axis
  xlab("Soft Threshold (power)") +  # These next three lines set the labels for our axes, and makes the formatting of the graph nicer
  ylab("Scale Free Topology Model Fit, signed R^2") + 
  ggtitle("Scale Independence") + theme_classic()

# We want to choose a value for our "power" parameter that is above the red line on this plot.
# BUT we don't want our power parameter to be excessively large. It's best to choose one at an inflection point, where increasing power further only gives minimal increases in signed R^2.
# BEFORE looking at the next block of code, what do you think our power parameter should be? 

# Once we've identified a power parameter, we can run blockwiseModules(), which creates our gene co-expression modules
# It groups genes based on similarity of expression patterns across all samples.
# This code might take a little while to run
bwnet <- blockwiseModules(normalized_counts, # This line passes our "normalized_counts" dataset that we created on lines 191-192 to the blockwiseModules function
                          maxBlockSize = 5000, # We will split the data into chunks of 5000 genes because our computer can't do all 20,000+ at once
                          TOMType = "signed", # Type of matrix we have
                          power = 14, # I chose a power parameter of 14, based on the plot we generated in the previous step. Do you understand why?
                          numericLabels = T, # Modules will be numbered
                          randomSeed = 1234) # Since this process is a random search, we won't necessarily get the same results every time. If we set a seed though, these results will be reproducible. 

# Write these results to a file in case we want to return to them
readr::write_rds(bwnet, file = "WGCNA_results.RDS")  

# In our "bwnet" object (the results of the WGCNA), there is a lot of information
# What we are most interested in is the "eigengene" for each sample in each module
# An eigengene is a "summary" value of the expression of all genes that are included in a single module (similar to an eigenvalue generated during a PCA)
# For example, if module 1 has 100 genes, we can summarize the expression of all 100 genes in the module into a single "eigengene" value of expression
# Each sample (Apoly_Dec_1, Apoly_Dec_2, etc.) will have a different eigengene value for a particular module, because each sample has a slightly different expression pattern
# For more info on eigengenes, take a look at section 4.9 of the tutorial linked at the top of this section
# Extract the eigengene values for each module into an object called module_eigengenes
module_eigengenes <- bwnet$MEs
# Take a look at the eigengene values for each sample and each module
head(module_eigengenes)

# Which modules have the biggest differences across our treatment groups?
# I.e., which clusters of co-expressed genes are most variable between December and February samples?

# First, let's make sure our samples are still in the correct order
all.equal(coldata$fishID, rownames(module_eigengenes))

# We'll create a "design" matrix based on the "treatment" variable so the program knows that's what we want to split our samples by
des_mat <- model.matrix(~ coldata$treatment)

# Now we're going to create a linear model on each module of co-expressed genes, to determine how the eigengene expression value relates to the treatment
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
fit <- limma::eBayes(fit)

# Apply a correction for multiple tests (we're testing a bunch of different modules), then obtain corrected stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

# Take a look at the results of our statistical test - modules that are most significantly related to "treatment" (i.e., show biggest differences in expression between Dec and Feb)
# will be at the top of the table
head(stats_df)

# Which module seems to be the MOST differentially expressed across our treatment groups?
# Let's investigate this module.

# First, let's just visualize what the eigengene for this module looks like between different treatment groups
# Set up the module eigengene using the sample labels that we need
module_1_df <- module_eigengenes %>%
  tibble::rownames_to_column("FishName") %>%
  dplyr::inner_join(coldata %>%
                      dplyr::select(fishID, treatment),
                    by = c("FishName" = "fishID"))


# Now plot the module
ggplot(module_1_df, aes(x=treatment, y=ME1, color=treatment)) + 
  geom_boxplot(width = 0.2, outlier.shape = NA) + 
  ggforce::geom_sina(maxwidth = 0.3) + theme_classic()

# In this case, we've got the treatment group on the x-axis, and the module eigengene values for each sample on the y-axis.
# What conclusion can we draw about this module of co-expressed genes between our treatment groups?

# Let's find out what genes are a part of module 1! Maybe we can understand the functionality of this co-expressed module!
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))
gene_module_key %>%
  dplyr::filter(module == "ME1")

# If each row in this table is a gene in module 1, how many genes are there in this module of co-expression?

# Let's write these results to a table so we can examine them more easily.
readr::write_tsv(gene_module_key, file = "WGCNA_module_1.tsv")
# Search through this results file like how we looked at the full output of DESeq. Are there any interesting genes in this module of co-expression?

# Finally, let's make a heatmap that summarizes the expression patterns of each gene in this module, across all samples. 
# First, we need to define a custom "function" (or command) in R that will make building a heatmap straightforward. 
# I got this code straight from the tutorial linked at the top of this section. 
make_module_heatmap <- function(module_name, expression_mat = normalized_counts, # These next few lines specify the function name (make_module_heatmap) and any arguments the function takes
                                metadata_df = coldata, gene_module_key_df = gene_module_key, # as well as any default values of those arguments. For example, it requires a module_name, as well as an expression matrix
                                module_eigengenes_df = module_eigengenes) { # There's no default module name, but the default expression matrix name is normalized_counts, etc.
  
  module_eigengene <- module_eigengenes_df %>% # This block of code creates a little table with the fish ID and the module eigengene value associated with each fish
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("fishID")
  
  col_annot_df <- metadata_df %>% # This block of code adds the treatment information (i.e., Dec or Feb) to that table
    dplyr::select(treatment, fishID) %>%
    dplyr::inner_join(module_eigengene, by = "fishID") %>%
    dplyr::arrange(treatment, fishID) %>%
    tibble::column_to_rownames("fishID")
  
  col_annot <- ComplexHeatmap::HeatmapAnnotation( # This block of code creates a small barplot below the actual heatmap, which will contain the eigengene (i.e., summary) expression values for each sample. The samples will be blocked and colored by treatment.
    treatment = col_annot_df$treatment,
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    col = list(treatment = c("December" = "blue2", "February" = "coral2"))
  )
  
  module_genes <- gene_module_key_df %>% # This creates an object called module_genes, which selects all of the gene names in our module of interest
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  mod_mat <- expression_mat %>% # This grabs the expression values for each gene in our module
    t() %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% module_genes) %>%
    dplyr::select(rownames(col_annot_df)) %>%
    as.matrix()
  
  mod_mat <- mod_mat %>% # This normalizes the expression values
    t() %>%
    scale() %>%
    t()
  
  color_func <- circlize::colorRamp2( # And this creates our colorization scale
    c(-2, 0, 2),
    c("blue2", "white", "coral2")
  )
  
  heatmap <- ComplexHeatmap::Heatmap(mod_mat, # Finally, this plots the expression values of each gene in each sample on a heatmap
                                     name = module_name,
                                     col = color_func,
                                     bottom_annotation = col_annot,
                                     cluster_columns = F,
                                     show_row_names = F,
                                     show_column_names = F)
  
  return(heatmap) # And the function will output the heatmap
  
}


module_1_heatmap <- make_module_heatmap(module_name = "ME1")
module_1_heatmap

# What does this heatmap tell us??


#################     Second sample dataset      #################

# Use the code from above, as well as what we discussed in class today, to analyze another sample dataset. 
# The files you'll need for this analysis are on our course GitHub page.
# The read count data is in the file: count_data_2.csv, and the sample information is in: col_data_2.csv
# You'll notice that none of the genes in this dataset have full names, they're all just called TRINITY_DNXXXX.
# You'll have to search these TRINITY ID numbers in the file annotation_info.csv (also on the GitHub page) to figure out what genes are what, when you get your results.

# This second sample dataset is from the example that I discussed earlier in class today, about pinfish thermal stress.
# Briefly, we had three treatment groups: control (sample names begin with C), moderate heat stress (sample names begin with M), and extreme heat stress (sample names begin with H)
# The samples included in this dataset contain liver tissue from each fish. 
# The information in col_data_2.csv tells you which treatment group a sample was in (21 - control, 24 - moderate, or 27 - extreme)

