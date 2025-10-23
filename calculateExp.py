#!/usr/bin/env python3
"""
eModule: Enhancer Analysis Framework - Expression Quantification Module

This script quantifies enhancer expression levels from RNA-seq bigWig files by 
counting reads mapped within each enhancer region. It serves as the first step 
in the eModule pipeline for enhancer-mediated gene regulation analysis.
"""

import os
import sys, getopt
import pyBigWig


def read_enhancer(annoFile):
    """
    Read enhancer annotation file and extract genomic coordinates.
    
    Parameters:
    -----------
    annoFile : str
        Path to enhancer annotation file in BED format (chromosome\tstart\tend)
    
    Returns:
    --------
    list of tuples
        List containing (chromosome, start, end, locus) for each enhancer
        where locus is formatted as "chr:start-end"
    """
    enhancers = []
    with open(annoFile) as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            locus = '%s:%s-%s' % (chrom, start, end)
            start = int(start)
            end = int(end)
            enhancers.append((chrom,
                              start,
                              end,
                              locus))
    return enhancers
    
    
def extract_exp(bwfile, enhancers, outFile):
    """
    Extract enhancer expression values from bigWig file.
    
    This function calculates the average read coverage across each enhancer region
    using pyBigWig library. The coverage values represent the transcriptional 
    activity of enhancer RNAs (eRNAs).
    
    Parameters:
    -----------
    bwfile : str
        Path to bigWig file containing RNA-seq coverage data
    enhancers : list
        List of enhancer tuples from read_enhancer()
    outFile : str
        Output file path to write expression values
    
    Output format:
    -------------
    locus\taverage_coverage
    chr1:1000-2000\t15.67
    
    Notes:
    ------
    - Uses exact=True for precise coverage calculation
    - Results are rounded to 2 decimal places for consistency
    - Each enhancer's expression is quantified as mean read coverage per base
    """
    # Open bigWig file for reading
    bw = pyBigWig.open(bwfile)
    counts = {}
    
    # Calculate average coverage for each enhancer
    for chrom, start, end, locus in enhancers:
    # Get average value over the enhancer range
        count = bw.stats(chrom, start, end, exact=True)[0]
        counts[locus] = str(round(count, 2))
    bw.close()
	
    # Write results to output file
    with open(outFile, 'w') as f_re:
        for t in enhancers:
            locus = t[-1]
            count = counts[locus]
            f_re.write('\t'.join([locus, count])+'\n')


def usage():
    """
    Print usage information and parameter descriptions.
    """
    print("""Parameters:
        --annoFile       enhancer_annotation_file. 
        --bwfile         bigwig file. This file can be downloaded from the Recount3 platform.
        --outFile        the file to write result.
        """)
    print()
    print('Example: python calculateExp.py  --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene --bwfile gtex.base_sums.ADIPOSE_TISSUE_GTEX-1A3MV-2126-SM-718BV.1.ALL.bw --outFile rawExp')
    
       
if __name__ == '__main__':
    # Parse command line arguments
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['annoFile=', 'bwfile=', 'outFile='])
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
        elif opt == '--bwfile':
            bwfile = str(arg)
        elif opt == '--outFile':
            outFile = str(arg)
    
    # Main workflow
    print("Reading enhancer annotations...")
    enhancers = read_enhancer(annoFile)
    print(f"Found {len(enhancers)} enhancers")
    
    print("Extracting expression values...")
    extract_exp(bwfile, enhancers, outFile)
    print("Expression quantification completed!")