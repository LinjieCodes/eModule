#!/usr/bin/env python3
"""
eModule: Enhancer Analysis Framework - Gene Regulation Module Identification

This script implements the core algorithm for identifying enhancer-mediated 
gene regulatory modules through co-expression analysis and mediation testing.
It integrates enhancer expression, gene expression, and transcription factor 
binding data to infer regulatory relationships.

The algorithm follows these key steps:
1. Identifies co-expressed triplets (TF-enhancer-gene) using Spearman correlation
2. Filters for enhancer-mediated regulation by ensuring TF binds enhancer 
   but not gene promoter
3. Performs bootstrap mediation analysis to test indirect effects
4. Identifies significant enhancer-mediated regulatory modules
"""

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
import numpy as np
from scipy import stats
import pandas as pd
import re
import sys, getopt


def read_enhancer_exp(eExpCsv):
    """
    Read enhancer expression matrix from CSV file.
    
    Parameters:
    -----------
    eExpCsv : str
        Path to enhancer expression CSV file (enhancers × samples)
    
    Returns:
    --------
    pd.DataFrame
        Enhancer expression matrix with enhancers as rows, samples as columns
    """
    eExp_df = pd.read_csv(eExpCsv, index_col=0)
    return eExp_df
    
    
def read_gene_exp(gExpCsv):
    """
    Read gene expression matrix from CSV file.
    
    Parameters:
    -----------
    gExpCsv : str
        Path to gene expression CSV file (genes × samples)
    
    Returns:
    --------
    pd.DataFrame
        Gene expression matrix with genes as rows, samples as columns
    """
    gExp_df = pd.read_csv(gExpCsv, index_col=0)
    return gExp_df
    
    
def merge_exp(eExp_df, gExp_df):
    """
    Merge enhancer and gene expression matrices for co-expression analysis.
    
    Parameters:
    -----------
    eExp_df : pd.DataFrame
        Enhancer expression matrix
    gExp_df : pd.DataFrame
        Gene expression matrix
    
    Returns:
    --------
    pd.DataFrame
        Combined expression matrix with both enhancers and genes as rows
    """
    # Find common samples across both matrices
    common_columns = eExp_df.columns.intersection(gExp_df.columns)
    merged_exp_df = pd.concat([eExp_df[common_columns], gExp_df[common_columns]])
    
    # Log-transform to stabilize variance (add pseudocount of 1)
    merged_exp_df = np.log2(merged_exp_df + 1)
    return merged_exp_df
    
  
def read_eTFBS(eTFBS_file, tfcutoff):
    """
    Read transcription factor binding sites (TFBS) for enhancers.
    
    Parameters:
    -----------
    eTFBS_file : str
        Path to enhancer TFBS file (enhancer\tTFs\tscore format)
    tfcutoff : int
        Minimum binding score threshold for considering TFBS significant
    
    Returns:
    --------
    dict
        Mapping of enhancer loci to sets of bound transcription factors
    """
    eTFBS = {}
    with open(eTFBS_file) as f:
        f.readline()    # Skip header
        for line in f:
            enhancer, tfs, score = line.strip().split('\t')
            if int(score) > tfcutoff:
                for tf in tfs.split('::'):
                    if enhancer not in eTFBS:
                        eTFBS[enhancer] = set()
                    eTFBS[enhancer].add(tf.upper())
    return eTFBS
    
    
def read_gTFBS(gTFBS_file, tfcutoff):
    """
    Read transcription factor binding sites (TFBS) for gene promoters.
    
    Parameters:
    -----------
    gTFBS_file : str
        Path to gene TFBS file (gene\tTFs\tscore format)
    tfcutoff : int
        Minimum binding score threshold for considering TFBS significant
    
    Returns:
    --------
    dict
        Mapping of gene symbols to sets of bound transcription factors
    """
    gTFBS = {}
    with open(gTFBS_file) as f:
        f.readline()  # Skip header
        for line in f:
            gene, tfs, score = line.strip().split('\t')
            if int(score) > tfcutoff:
                for tf in tfs.split('::'):
                    if gene not in gTFBS:
                        gTFBS[gene] = set()
                    gTFBS[gene].add(tf.upper())
    return gTFBS
    

def read_enhancer(eFile):
    """
    Read list of enhancers to analyze.
    
    Parameters:
    -----------
    eFile : str
        Path to file containing enhancer loci (one per line)
    
    Returns:
    --------
    list
        List of enhancer loci to analyze
    """
    enhancers = []
    with open(eFile) as f:
        for line in f:
            enhancers.append(line.strip())
    return enhancers

    
def obtain_near_gene(enhancers, gtfFile):
    """
    Identify genes within 1Mb of each enhancer using GTF annotation.
    
    This function implements a genomic proximity search to find candidate 
    target genes for each enhancer based on linear distance.
    
    Parameters:
    -----------
    enhancers : list
        List of enhancer loci to analyze
    gtfFile : str
        Path to GTF annotation file containing gene coordinates
    
    Returns:
    --------
    dict
        Mapping of enhancer loci to sets of nearby gene symbols (within 1Mb)
    """
    genes = {}
    
    # Parse GTF file to extract protein-coding genes
    with open(gtfFile) as f:
        for line in f:
            if '\tgene\t' in line and 'gene_biotype "protein_coding"' in line:
                cols = line.split('\t')
                chrom = 'chr'+cols[0]
                start = int(cols[3])
                region_start = start - 1000000  # 1Mb upstream
                region_end = start + 1000000    # 1Mb downstream
                
                # Extract gene symbol
                gene_symbol = re.search('gene_name "(.+?)"', line).group(1)
                if chrom not in genes:
                    genes[chrom] = []
                genes[chrom].append([region_start, region_end, gene_symbol])
    
    # Sort genes on each chromosome by position    
    for chrom in genes:
        genes[chrom].sort()
    
    # Find genes within 1Mb of each enhancer    
    near_genes = {}
    for enhancer in enhancers:
        chrom, locus = enhancer.split(':')
        eStart = int(locus.split('-')[0])
        for region_start, region_end, gene_symbol in genes[chrom]:
            # Skip if gene region is before enhancer
            if region_end < eStart:
                continue
            # Stop if gene region is after enhancer
            elif region_start > eStart:
                break
            # Gene is within 1Mb of enhancer
            elif region_start < eStart < region_end:
                if enhancer not in near_genes:
                    near_genes[enhancer] = set()
                near_genes[enhancer].add(gene_symbol)
    return near_genes
 
 
def boot_mediation(tf, enh, gene, X, M, Y,
                   n_boot=10000, seed=42, n_jobs=-1):
    """
    Perform bootstrap mediation analysis to test TF → enhancer → gene regulation.
    
    This function implements non-parametric bootstrap mediation testing to assess
    whether the relationship between TF and gene is mediated by enhancer expression.
    
    Parameters:
    -----------
    tf : str
        Transcription factor symbol
    enh : str  
        Enhancer locus
    gene : str
        Target gene symbol
    X : array-like
        TF expression values
    M : array-like
        Enhancer expression values  
    Y : array-like
        Gene expression values
    n_boot : int
        Number of bootstrap resamples (default: 10000)
    seed : int
        Random seed for reproducibility
    n_jobs : int
        Number of parallel jobs for bootstrap sampling
    
    Returns:
    --------
    dict
        Dictionary containing mediation analysis results:
        - TF, Enhancer, Target_gene: Identifiers
        - ACME: Average causal mediation effect (indirect effect)
        - CI_lower, CI_upper: 95% bootstrap confidence interval
        - Proportion Mediated: Indirect effect / total effect
        - P_val: Bootstrap p-value for mediation significance
    """
    # Original estimate
    # M ~ X (path a)
    a_hat = sm.OLS(M, sm.add_constant(X)).fit().params[1]
    
    # Y ~ X + M (paths b and c')
    full = sm.OLS(Y, sm.add_constant(pd.DataFrame({'X': X, 'M': M}))).fit()
    b_hat = full.params['M']  # path b
    c_hat = full.params['X']  # path c'
    
    # Indirect effect (a * b)
    indirect_hat = a_hat * b_hat
    mediated_prop = indirect_hat/c_hat

    # Bootstrap sampling
    rng = np.random.RandomState(seed)
    
    def _iter(i):
        # Resample with replacement
        idx = rng.choice(len(X), size=len(X), replace=True)
        Xb, Mb, Yb = X[idx], M[idx], Y[idx]
        
        # Calculate indirect effect for bootstrap sample
        ab = (sm.OLS(Mb, sm.add_constant(Xb)).fit().params[1] *
              sm.OLS(Yb, sm.add_constant(pd.DataFrame({'X': Xb, 'M': Mb}))).fit().params['M'])
        return ab

    # Parallel bootstrap sampling
    boot_ests = np.array(Parallel(n_jobs=n_jobs)(
        delayed(_iter)(i) for i in range(n_boot)))

    # Confidence interval and p-value (percentile bootstrap)
    ci_lower, ci_upper = np.percentile(boot_ests, [2.5, 97.5])
    p_boot = 2 * min(np.mean(boot_ests <= 0),
                     np.mean(boot_ests >= 0))

    return {'TF': tf, 'Enhancer': enh, 'Target_gene': gene,
            'ACME(average causal mediation effect)': indirect_hat,
            'CI_lower': ci_lower, 'CI_upper': ci_upper,
            'Proportion Mediated(ACME/Total Effect)':mediated_prop,
            'P_val': p_boot}

            
def identify_targets(enhancers,
                     near_genes,
                     merged_exp_df,
                     eTFBS,
                     gTFBS,
                     corr_cutoff,
                     pval_cutoff,
                     n_boot,
                     outFile):
    """
    Identify enhancer-mediated gene regulatory modules.
    
    This function implements the core algorithm:
    1. For each enhancer, identify regulating TFs (bind enhancer, co-expressed)
    2. For each regulating TF, find candidate genes within 1Mb
    3. Test co-expression between TF-enhancer-gene triplets
    4. Filter for enhancer-mediated regulation (TF binds enhancer but not gene promoter)
    5. Perform mediation analysis to confirm indirect effects
    
    Parameters:
    -----------
    enhancers : list
        List of enhancer loci to analyze
    near_genes : dict
        Mapping of enhancers to nearby genes
    merged_exp_df : pd.DataFrame
        Combined enhancer and gene expression matrix
    eTFBS : dict
        Enhancer-TF binding information
    gTFBS : dict
        Gene-TF binding information
    corr_cutoff : float
        Spearman correlation coefficient threshold
    pval_cutoff : float
        P-value threshold for significance testing
    outFile : str
        Output file path for regulatory modules
    """
    triplets_dict = {}
    
    print("Identifying co-expressed TF-enhancer-gene triplets...")
    for enhancer in enhancers:
        if enhancer not in eTFBS or enhancer not in near_genes:
            continue
            
        # Find regulating TFs that bind this enhancer and are co-expressed
        regulating_TFs = set()
        for eTF in eTFBS[enhancer]:
            if eTF in merged_exp_df.index:
                r, pval = stats.spearmanr(merged_exp_df.loc[enhancer],
                                          merged_exp_df.loc[eTF])
                if abs(r) > corr_cutoff and pval < pval_cutoff:
                    regulating_TFs.add(eTF)
        
        # Find target genes for each regulating TF
        target_genes = {}
        if regulating_TFs:
            for nearGene in near_genes[enhancer]:
                for regulating_TF in regulating_TFs:
                    # Check if TF does NOT bind gene promoter (enhancer-mediated regulation)
                    if nearGene not in gTFBS or regulating_TF not in gTFBS[nearGene]:
                        if nearGene in merged_exp_df.index:
                            # Test co-expression between all three components
                            r1, pval1 = stats.spearmanr(merged_exp_df.loc[enhancer],
                                                        merged_exp_df.loc[nearGene])
                            r2, pval2 = stats.spearmanr(merged_exp_df.loc[regulating_TF],
                                                        merged_exp_df.loc[nearGene])
                            
                            if (abs(r1) > corr_cutoff and pval1 < pval_cutoff and 
                            abs(r2) > corr_cutoff and pval2 < pval_cutoff):
                                if regulating_TF not in target_genes:
                                    target_genes[regulating_TF] = set()
                                target_genes[regulating_TF].add(nearGene)                                
        
        if target_genes:
            triplets_dict[enhancer] = target_genes
    
    # Convert triplets to list for analysis
    triplets_list = []
    for enhancer in triplets_dict:
        for TF in triplets_dict[enhancer]:
            for gene in triplets_dict[enhancer][TF]:
               triplets_list.append((TF, enhancer, gene))
               
    print(f"Testing {len(triplets_list)} triplets with mediation analysis...")
    
    # Perform bootstrap mediation analysis in parallel
    results = Parallel(n_jobs=-1)(
        delayed(boot_mediation)(TF, enhancer, gene,
                                merged_exp_df.loc[TF], 
                                merged_exp_df.loc[enhancer], 
                                merged_exp_df.loc[gene])
        for (TF, enhancer, gene) in triplets_list
    )

    # Process mediation results
    df_boot = pd.DataFrame(results)
    df_boot['FDR'] = multipletests(df_boot['P_val'], method='fdr_bh')[1]

    # Filter for significant mediation effects
    sig = df_boot[(df_boot['FDR'] < pval_cutoff) &
                  (df_boot['CI_lower'] * df_boot['CI_upper'] > 0)]
    
    print(f"Found {len(sig)} significant enhancer-mediated regulatory relationships")
    
    # Organize results into regulatory modules
    modules = {}
    for idx, row in sig.iterrows():
        tf = row['TF']
        enhancer = row['Enhancer']
        gene = row['Target_gene']
        
        if enhancer not in modules:
            modules[enhancer] = {}
        if tf not in modules[enhancer]:
            modules[enhancer][tf] = set()
        modules[enhancer][tf].add(gene)
    
    # Write regulatory modules to output file    
    with open(outFile, 'w') as f:
        f.write('Enhancer\tRegulating TF\tTarget genes\n')
        for enhancer in modules:
            for tf in modules[enhancer]:
                target_str = ', '.join(modules[enhancer][tf])
                f.write('\t'.join([enhancer, tf, target_str])+'\n')
                

def usage():
    """
    Print usage information and parameter descriptions.
    """
    print("""Parameters:
        --eExpCsv         enhancer expression file (in csv format). 
        --gExpCsv         gene expression file (in csv format). 
        --eTFBS_file      Enhancers' TFBS profiles.
        --gTFBS_file      Genes' TFBS profiles.
        --eFile           Enhancers to analyze.
        --gtfFile         Gene GTF file. This file can be downloaded from Ensembl.
        --n_boot          specifies the number of resampling times.
        --tfcutoff        TF binding score cutoff.
        --rcutoff         Correlation coefficient cutoff.
        --pcutoff         P-value cutoff.
        --outFile         Output filename.
        """)
    print()
    print('Example: python identifyModule.py  --eExpCsv Spleen_RPM.csv  --gExpCsv Spleen_geneExp.csv --eTFBS_file Enhancer_TFBS_demo --gTFBS_file Promoter_TFBS_demo --eFile SpleenSexBiasedEnhancer --gtfFile Homo_sapiens.GRCh38.101.gtf --n_boot 10000  --tfcutoff 400 --rcutoff 0.3 --pcutoff 0.05  --outFile module_out')
    
            
if __name__ == '__main__':
    # Parse command line arguments
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['eExpCsv=',
                                                 'gExpCsv=', 
                                                 'eTFBS_file=', 
                                                 'gTFBS_file=', 
                                                 'eFile=',
                                                 'n_boot=',                                                 
                                                 'gtfFile=', 
                                                 'tfcutoff=',
                                                 'rcutoff=', 
                                                 'pcutoff=', 
                                                 'outFile='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # Process command line options with defaults
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '--eExpCsv':
            eExpCsv = str(arg)
        elif opt == '--gExpCsv':
            gExpCsv = str(arg)
        elif opt == '--eTFBS_file':
            eTFBS_file = str(arg)
        elif opt == '--gTFBS_file':
            gTFBS_file = str(arg)
        elif opt == '--eFile':
            eFile = str(arg)
        elif opt == '--gtfFile':
            gtfFile = str(arg)
        elif opt == '--n_boot':
            n_boot = str(arg)
        elif opt == '--rcutoff':
            corr_cutoff = float(arg)
        elif opt == '--tfcutoff':
            tfcutoff = int(arg)
        elif opt == '--pcutoff':
            pval_cutoff = float(arg)
        elif opt == '--outFile':
            outFile = str(arg)
    
    # Set default values if not provided
    if 'corr_cutoff' not in locals():
        corr_cutoff = 0.3
    if 'tfcutoff' not in locals():
        tfcutoff = 400
    if 'pval_cutoff' not in locals():
        pval_cutoff = 0.05
    if 'n_boot' not in locals():
        n_boot = 10000
    
    print("Starting enhancer module identification...")
    print(f"Parameters: corr_cutoff={corr_cutoff}, tfcutoff={tfcutoff}, pval_cutoff={pval_cutoff}")
    
    # Load data
    print("Loading expression data...")
    eExp_df = read_enhancer_exp(eExpCsv)
    gExp_df = read_gene_exp(gExpCsv)
    
    print("Merging expression matrices...")
    merged_exp_df = merge_exp(eExp_df, gExp_df)
    print(f"Combined expression matrix: {merged_exp_df.shape}")
    
    print("Loading TF binding data...")
    eTFBS = read_eTFBS(eTFBS_file, tfcutoff)
    gTFBS = read_gTFBS(gTFBS_file, tfcutoff)
    
    print("Reading enhancer list...")
    enhancers = read_enhancer(eFile) 
    print(f"Analyzing {len(enhancers)} enhancers")
    
    print("Finding nearby genes...")
    near_genes = obtain_near_gene(enhancers, gtfFile)
    print(f"Found nearby genes for {len(near_genes)} enhancers")
    
    # Run module identification
    identify_targets(enhancers,
                     near_genes,
                     merged_exp_df,
                     eTFBS,
                     gTFBS,
                     corr_cutoff,
                     pval_cutoff,
                     n_boot,
                     outFile)
    
    print("Module identification completed!")