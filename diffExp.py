#!/usr/bin/env python3
"""
eModule: Enhancer Analysis Framework - Differential Expression Analysis Module

This script identifies differentially expressed enhancers using a generalized 
linear model (GLM) approach. It specifically focuses on sex-biased enhancer 
expression while controlling for multiple confounding factors including 
demographic variables, technical covariates, and population structure.

The script implements a robust statistical framework that:
1. Controls for gene expression PCs that correlate with sex
2. Uses Poisson GLM appropriate for count data
3. Applies multiple testing correction (Benjamini-Hochberg FDR)
4. Identifies enhancers with significant sex-biased expression
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import glm
from statsmodels.genmod.families import Poisson
from statsmodels.genmod.families.links import identity
from scipy.stats import pointbiserialr
import statsmodels
import os
import sys, getopt


def diff_exp(expFile,
             sampleFile,
             tissue,
             outFile):
    """
    Perform differential expression analysis to identify sex-biased enhancers.
    
    This function implements a comprehensive statistical pipeline:
    - Controls for sex-correlated gene expression PCs to avoid confounding
    - Fits Poisson GLM models with multiple covariates
    - Calculates fold changes between sexes
    - Applies multiple testing correction
    - Filters for biologically significant differences
    
    Parameters:
    -----------
    expFile : str
        Path to enhancer expression matrix (CSV format, enhancers × samples)
    sampleFile : str
        Path to sample attribute file containing covariates
    tissue : str
        Tissue identifier for output labeling
    outFile : str
        Output file path for results
    
    Statistical Model:
    ------------------
    Uses Poisson GLM with formula:
    expression ~ Sex + Age + GenotypePCs + GeneExpressionPCs + RIN + PMI
    
    Where Sex is coded as 1=male, 2=female
    
    Output Columns:
    ---------------
    Tissue: Tissue name
    Enhancer: Enhancer locus (chr:start-end)
    Male_exp: Median log2-CPM expression in males
    Female_exp: Median log2-CPM expression in females  
    log2(FC): Log2 fold change (male - female)
    coef: GLM coefficient for sex effect
    Pval: Raw p-value from GLM
    FDR: Benjamini-Hochberg adjusted p-value
    """
    results = []
    p_values = []

    # Read expression data (enhancers as columns, samples as index)
    print("Reading expression data...")
    expression_data = pd.read_csv(expFile, index_col=0)  # Transpose to samples × enhancers
    expression_data = expression_data.T
    print(f"Expression matrix shape: {expression_data.shape}")
    
    # Read sample attributes
    print("Reading sample attributes...")
    individual_attributes = pd.read_csv(sampleFile, index_col=0)
    individual_attributes = individual_attributes.T
    print(f"Sample attribute matrix shape: {individual_attributes.shape}")
    
    # Identify and exclude gene expression PCs that correlate with sex
    # This prevents removal of sex-associated biological signals
    print("Identifying sex-correlated expression PCs...")
    correlatedPCs = set()
    expPCs = [PC for PC in individual_attributes.columns if 'InferredCov' in PC]
    
    for expPC in expPCs:
        corr, p_value = pointbiserialr(individual_attributes['Sex'], individual_attributes[expPC])
        if p_value < 0.05:
            correlatedPCs.add(expPC)

    print(f"Excluding {len(correlatedPCs)} sex-correlated PCs")
    
    # Remove sex-correlated gene expression PCs to avoid over-correction
    for expPC in correlatedPCs:
        individual_attributes = individual_attributes.drop(expPC, axis = 1)
    
    # Construct GLM formula with all available covariates
    formula = 'expression ~ ' + ' + '.join(individual_attributes.columns)
    
    # Process each enhancer
    print("Fitting GLM models for each enhancer...")
    for enhancer in expression_data.columns:
        # Get expression data for current enhancer
        enhancer_expression_data = expression_data[[enhancer]]
        enhancer_expression_data.columns  = ['expression']
        
        # Calculate median expression by sex (with pseudocount)
        male_median = round(float(np.median(np.log2(enhancer_expression_data[individual_attributes['Sex']==1]+0.01))), 3)
        female_median = round(float(np.median(np.log2(enhancer_expression_data[individual_attributes['Sex']==2]+0.01))), 3)
        
        # Calculate log2 fold change (male - female)
        fold_change = round(male_median - female_median, 3)
    
        # Combine expression and covariate data
        combined_data = enhancer_expression_data.join(individual_attributes)
    
        # Fit Poisson GLM
        model = glm(formula, combined_data, family=Poisson()).fit()

        # Extract sex coefficient and p-value
        model_summary = model.summary2().tables[1]
        pval = model_summary['P>|z|']['Sex']
        coef = round(model_summary['Coef.']['Sex'], 3)
        
        p_values.append(pval)
        results.append([tissue, enhancer, male_median, female_median, fold_change, coef, format(pval, '.3e')])

    # Apply multiple testing correction
    print("Applying Benjamini-Hochberg correction...")
    p_values = np.array(p_values)
    r, fdrs, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(p_values, alpha=0.05, method='fdr_bh')
    
    # Add FDR values to results
    for i in range(len(fdrs)):
        results[i].append(fdrs[i])

    # Write significant results to output file
    print("Writing results...")                
    with open(outFile, 'w') as f:
        # Write header
        f.write('\t'.join(['Tissue',
                           'Enhancer',
                           'Male_exp',
                           'Female_exp',
                           'log2(FC)',
                           'coef',
                           'Pval',
                           'FDR'])+'\n')
        
        # Write significant enhancers
        significant_count = 0
        for tissue, enhancer, male_median, female_median, fold_change, coef, pval, fdr in results:
            # Filter for significant sex-biased enhancers
            # FDR < 0.05, |log2FC| > 1, and consistent direction between coef and FC
            if fdr < 0.05 and abs(fold_change)>1 and coef*fold_change<0:
                f.write('\t'.join([tissue,
                                   enhancer,
                                   str(male_median),
                                   str(female_median),
                                   str(fold_change),
                                   str(coef),
                                   str(pval),
                                   str(format(fdr, '.3e'))])+'\n')
                significant_count += 1 

    print(f"Analysis complete! Found {significant_count} sex-biased enhancers")                


def usage():
    """
    Print usage information and parameter descriptions.
    """
    print("""Parameters:
        --expFile        enhancer RPM matrix. 
        --sampleFile     sample_attribute_file of the tissue. The file can be downloaded from Recount3 platform.  
        --tissue         the tissue label.
        --outFile        the file to write result.
        """)
    print()
    print('Example: python diffExp.py  --expFile Spleen_RPM.csv --sampleFile Spleen_sample.csv --tissue Spleen --outFile sexBiasedEnhancer')
    

if __name__ == '__main__':
    # Parse command line arguments
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['expFile=', 'sampleFile=', 'tissue=', 'outFile='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    # Process command line options
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '--expFile':
            expFile = str(arg)
        elif opt == '--sampleFile':
            sampleFile = str(arg)
        elif opt == '--tissue':
            tissue = str(arg)
        elif opt == '--outFile':
            outFile = str(arg)  
    
    # Run differential expression analysis
    print(f"Starting differential expression analysis for {tissue}...")
    diff_exp(expFile,
             sampleFile,
             tissue,
             outFile)
    print("Analysis completed!")
