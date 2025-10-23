#!/usr/bin/env Rscript
"""
eModule: Enhancer Analysis Framework - TMM Normalization R Script

This R script provides TMM (Trimmed Mean of M-values) normalization using 
the edgeR package. TMM normalization is a robust method for RNA-seq data 
that accounts for library size differences and composition bias.

The script is called from Python via rpy2 to perform normalization on 
enhancer expression count matrices.
"""

# Load required library
library(edgeR)

"""
Perform TMM normalization on count matrix

Parameters:
-----------
count : matrix or data.frame
    Raw count matrix with genes/enhancers as rows and samples as columns
output : character
    Output filename for normalized counts

Returns:
--------
None (writes normalized counts to file)

Notes:
------
- Uses edgeR's DGEList object for count data storage
- Applies calcNormFactors for TMM normalization
- Outputs normalized counts in CSV format
"""

tmm_nor <- function(count, output){
	# Create DGEList object from count matrix
	y <- DGEList(counts=count)
	
	# Calculate TMM normalization factors
	y <- calcNormFactors(y)
	
	# Write to file
	write.table(y$counts, sep =",", file=output)
}