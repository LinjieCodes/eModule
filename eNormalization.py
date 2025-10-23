#!/usr/bin/env python3
"""
eModule: Enhancer Analysis Framework - Expression Normalization Module

This script performs comprehensive normalization of enhancer expression data using 
TMM (Trimmed Mean of M-values) and CPM (Counts Per Million) normalization methods.
It processes raw enhancer expression counts across multiple tissues and samples,
producing normalized expression matrices suitable for downstream differential analysis.

The script integrates R's edgeR package via rpy2 for TMM normalization, followed by 
CPM transformation and tissue-specific filtering to ensure robust enhancer expression.
"""

import os
import sys, getopt
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri

# Source R script containing TMM normalization function
r.source('tmm.r')
pandas2ri.activate()


def read_enhancer(annoFile):
    """
    Read enhancer annotation file and create index mapping.
    
    Parameters:
    -----------
    annoFile : str
        Path to enhancer annotation file in BED format
    
    Returns:
    --------
    tuple
        enhancer_index: list of enhancer loci (chr:start-end)
        lengths: dict mapping enhancer loci to their lengths in base pairs
    """
    enhancer_index = []
    lengths = {}
    with open(annoFile) as f:
        for line in f:
            chr, start, end = line.strip().split('\t')
            length = int(end) - int(start) + 1
            enhancer = '%s:%s-%s' % (chr, start, end)
            enhancer_index.append(enhancer)
            lengths[enhancer] = length
    return enhancer_index, lengths


def read_sample(sampleFile):
    """
    Parse sample attribute file to organize samples by tissue type.
    
    Parameters:
    -----------
    sampleFile : str
        Path to sample attribute file with format: tissue\tsampleID\n
    Returns:
    --------
    tuple
        tissue_samples: dict mapping tissue names to lists of sample IDs
        all_samples: set containing all unique sample IDs
    """
    tissue_samples = {}
    all_samples = set()
    with open(sampleFile) as f:
        f.readline()  # Skip header
        for line in f:
            cols = line.split('\t')
            tissue = cols[0]
            sample = cols[1]
            if tissue not in tissue_samples:
                tissue_samples[tissue] = []
            tissue_samples[tissue].append(sample)
            all_samples.add(sample)
    return tissue_samples, all_samples
    

def read_enhancerExp(enhancer_index, 
                     lengths, 
                     all_samples, 
                     expFolder):
    """
    Read raw enhancer expression files and aggregate counts across samples.
    
    This function processes per-base read count files for each sample and tissue,
    converting them to total read counts per enhancer by multiplying coverage 
    by enhancer length.
    
    Parameters:
    -----------
    enhancer_index : list
        List of enhancer loci identifiers
    lengths : dict
        Mapping of enhancer loci to their lengths
    all_samples : set
        Set of all sample IDs to process
    expFolder : str
        Path to folder containing tissue subfolders with expression files
    
    Returns:
    --------
    dict
        Mapping of sample IDs to lists of read counts for each enhancer
    """
    sample_counts = {}
    for tissue in os.listdir(expFolder):
        print(f"Processing tissue: {tissue}")
        for sampleFile in os.listdir(expFolder+tissue):
            # Extract sample ID from filename (GTEX-XXXXX-XXXXX-SM-XXXXX)
            sampleID = sampleFile[sampleFile.find('GTEX-'):(sampleFile.find('.ALL.bw')-2)]
            if sampleID in all_samples:
                sample_counts[sampleID] = []
                this_sample_counts = {}
                
                # Read per-base coverage file
                with open(expFolder+tissue+'/'+sampleFile) as f:
                    for line in f:
                        cols = line.split('\t')
                        if cols[0] in lengths:
                            # Convert per-base coverage to total read counts
                            # Multiply coverage by enhancer length
                            this_sample_counts[cols[0]] = round(float(cols[1]) * lengths[cols[0]])
                
                for enhancer in enhancer_index:
                    sample_counts[sampleID].append(this_sample_counts[enhancer])
    return sample_counts


def normalize_exp(enhancer_index,
                  tissue_samples,
                  sample_counts,
                  outFolder):    
    """
    Perform comprehensive normalization of enhancer expression data.
    
    This function implements a multi-step normalization pipeline:
    1. Creates expression matrix (enhancers × samples)
    2. Applies TMM normalization using edgeR via rpy2
    3. Converts to CPM (Counts Per Million) for cross-sample comparability
    4. Filters lowly expressed enhancers (median CPM < 1)
    5. Exports tissue-specific expression matrices
    
    Parameters:
    -----------
    enhancer_index : list
        List of enhancer identifiers
    tissue_samples : dict
        Mapping of tissue names to sample lists
    sample_counts : dict
        Mapping of sample IDs to expression counts
    outFolder : str
        Output directory for normalized expression files
    """
    # Create DataFrame with enhancers as rows, samples as columns
    sampleList = list(sample_counts.keys())
    expValues = []
    for sample in sampleList:
        expValues.append(sample_counts[sample])
    
    # Create expression matrix
    counts_df = pd.DataFrame(expValues, columns=enhancer_index, index=sampleList)
    counts_df = counts_df.T

    # Convert to R dataframe for TMM normalization
    counts_df_r = pandas2ri.py2rpy(counts_df)

    # Apply TMM normalization using edgeR
    outFile = 'TMM_read_counts.csv'
    r.tmm_nor(counts_df_r, outFile)

    outFile = 'TMM_read_counts.csv'
    
    # Read TMM-normalized read counts
    counts_tmm = pd.read_csv(outFile, index_col=0)

    # Apply CPM normalization
    constant = 1000000
    counts_rpm = counts_tmm.apply(lambda x: round((x / x.sum()) * constant))
    counts_rpm = counts_rpm.astype(int)

    # Export tissue-specific expression matrices
    for tissue in tissue_samples:
        samples = tissue_samples[tissue] 
        counts_rpm_tissue = counts_rpm[samples]
        
        # Calculate median expression across samples
        medians = counts_rpm_tissue.median(axis=1)
        
        # Filter out lowly expressed enhancers (median CPM < 1)
        columns_to_drop = medians[medians < 1].index
        counts_rpm_tissue = counts_rpm_tissue.drop(index=columns_to_drop)
        
        # Write tissue-specific expression matrix
        counts_rpm_tissue.to_csv(outFolder + '%s.csv' % tissue, index=True)
        print(f"Processed {tissue}: {counts_rpm_tissue.shape[0]} enhancers retained")


def usage():
    """
    Print usage information and parameter descriptions.
    """
    print("""Parameters:
        --annoFile       enhancer_annotation_file. 
        --sampleFile     sample_attribute_file. 
        --expFolder      Path of enhancers' raw expression.
        --outFolder      Output path to write enhancers' normalized expression.
        """)
    print()
    print('Example: python eNormalization.py  --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene --sampleFile sampleAttributes --expFolder eRNA_perbase_average/ --outFolder eRNA_RPM/')

    
if __name__ == '__main__':
    # Parse command line arguments
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['annoFile=', 'sampleFile=', 'expFolder=', 'outFolder='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # Process command line options    
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '--annoFile':
            annoFile = str(arg)
        elif opt == '--sampleFile':
            sampleFile = str(arg)
        elif opt == '--expFolder':
            expFolder = str(arg)
        elif opt == '--outFolder':
            outFolder = str(arg)   
    
    # Create output directory if it doesn't exist
    if not os.path.exists(outFolder):
        os.mkdir(outFolder)
    
    # Main workflow
    print("Reading enhancer annotations...")
    enhancer_index, lengths = read_enhancer(annoFile)
    print(f"Loaded {len(enhancer_index)} enhancers")

    print("Reading sample information...")
    tissue_samples, all_samples = read_sample(sampleFile)
    print(f"Found {len(tissue_samples)} tissues with {len(all_samples)} total samples")
    
    print("Reading expression data...")
    sample_counts = read_enhancerExp(enhancer_index, lengths, all_samples, expFolder)
    print(f"Processed expression data for {len(sample_counts)} samples")

    print("Performing normalization...")
    normalize_exp(enhancer_index,
                  tissue_samples,
                  sample_counts,
                  outFolder)
    print("Normalization completed!")