from scipy import stats
import pandas as pd
import re
import sys, getopt


def read_enhancer_exp(eExpCsv):
    eExp_df = pd.read_csv(eExpCsv, index_col=0)
    return eExp_df
    
    
def read_gene_exp(gExpCsv):
    gExp_df = pd.read_csv(gExpCsv, index_col=0)
    return gExp_df
    
    
def merge_exp(eExp_df, gExp_df):
    common_columns = eExp_df.columns.intersection(gExp_df.columns)
    merged_exp_df = pd.concat([eExp_df[common_columns], gExp_df[common_columns]])
    return merged_exp_df
    
  
def read_eTFBS(eTFBS_file, tfcutoff):
    eTFBS = {}
    with open(eTFBS_file) as f:
        f.readline()
        for line in f:
            enhancer, tfs, score = line.strip().split('\t')
            if int(score) > tfcutoff:
                for tf in tfs.split('::'):
                    if enhancer not in eTFBS:
                        eTFBS[enhancer] = set()
                    eTFBS[enhancer].add(tf.upper())
    return eTFBS
    
    
def read_gTFBS(gTFBS_file, tfcutoff):
    gTFBS = {}
    with open(gTFBS_file) as f:
        f.readline()
        for line in f:
            gene, tfs, score = line.strip().split('\t')
            if int(score) > tfcutoff:
                for tf in tfs.split('::'):
                    if gene not in gTFBS:
                        gTFBS[gene] = set()
                    gTFBS[gene].add(tf.upper())
    return gTFBS
    

def read_enhancer(eFile):
    enhancers = []
    with open(eFile) as f:
        for line in f:
            enhancers.append(line.strip())
    return enhancers

    
def obtain_near_gene(enhancers, gtfFile):
    genes = {}
    with open(gtfFile) as f:
        for line in f:
            if '\tgene\t' in line and 'gene_biotype "protein_coding"' in line:
                cols = line.split('\t')
                chrom = 'chr'+cols[0]
                start = int(cols[3])
                region_start = start - 1000000
                region_end = start + 1000000
                gene_symbol = re.search('gene_name "(.+?)"', line).group(1)
                if chrom not in genes:
                    genes[chrom] = []
                genes[chrom].append([region_start, region_end, gene_symbol])
                
    for chrom in genes:
        genes[chrom].sort()
        
    near_genes = {}
    for enhancer in enhancers:
        chrom, locus = enhancer.split(':')
        eStart = int(locus.split('-')[0])
        for region_start, region_end, gene_symbol in genes[chrom]:
            if region_end < eStart:
                continue
            elif region_start > eStart:
                break
            elif region_start < eStart < region_end:
                if enhancer not in near_genes:
                    near_genes[enhancer] = set()
                near_genes[enhancer].add(gene_symbol)
    return near_genes
    
    
def identify_targets(enhancers,
                     near_genes,
                     merged_exp_df,
                     eTFBS,
                     gTFBS,
                     corr_cutoff,
                     pval_cutoff,
                     outFile):
    eModules = {}
    for enhancer in enhancers:
        regulating_TFs = set()
        if enhancer in eTFBS:
            for eTF in eTFBS[enhancer]:
                if eTF in merged_exp_df.index:
                    r, pval = stats.spearmanr(merged_exp_df.loc[enhancer],
                                              merged_exp_df.loc[eTF])
                    if abs(r) > corr_cutoff and pval < pval_cutoff:
                        regulating_TFs.add(eTF)
                    
        target_genes = {}
        if regulating_TFs:
            for nearGene in near_genes[enhancer]:
                for regulating_TF in regulating_TFs:
                    if nearGene not in gTFBS or regulating_TF not in gTFBS[nearGene]:
                        if nearGene in merged_exp_df.index:
                            r1, pval1 = stats.spearmanr(merged_exp_df.loc[enhancer],
                                                        merged_exp_df.loc[nearGene])
                            r2, pval2 = stats.spearmanr(merged_exp_df.loc[regulating_TF],
                                                        merged_exp_df.loc[nearGene])
                            if (abs(r1) > corr_cutoff and pval1 < pval_cutoff) and (abs(r2) > corr_cutoff and pval2 < pval_cutoff):
                                if regulating_TF not in target_genes:
                                    target_genes[regulating_TF] = set()
                                target_genes[regulating_TF].add(nearGene)                                
        
        if target_genes:
            eModules[enhancer] = target_genes
            
    with open(outFile, 'w') as f:
        f.write('Enhancer\tRegulating TF\tTarget genes\n')
        for enhancer in eModules:
            for regulating_TF in eModules[enhancer]:
                target_str = ', '.join(eModules[enhancer][regulating_TF])
                f.write('\t'.join([enhancer, regulating_TF, target_str])+'\n')
    

def usage():
    print("""Parameters:
        --eExpCsv         enhancer expression file (in csv format). 
        --gExpCsv         gene expression file (in csv format). 
        --eTFBS_file      Enhancers' TFBS profiles.
        --gTFBS_file      Genes' TFBS profiles.
        --eFile           Enhancers to analyze.
        --gtfFile         Gene GTF file. This file can be downloaded from Ensembl.
        --tfcutoff        TF binding score cutoff.
        --rcutoff         Correlation coefficient cutoff.
        --pcutoff         P-value cutoff.
        --outFile         Output filename.
        """)
    print()
    print('Example: python identifyModule.py  --eExpCsv  --gExpCsv --eTFBS_file --gTFBS_file --eFile --gtfFile --tfcutoff --rcutoff --pcutoff --outFile')
    
            
if __name__ == '__main__':
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['eExpCsv=',
                                                 'gExpCsv=', 
                                                 'eTFBS_file=', 
                                                 'gTFBS_file=', 
                                                 'eFile=', 
                                                 'gtfFile=', 
                                                 'tfcutoff=',
                                                 'rcutoff=', 
                                                 'pcutoff=', 
                                                 'outFile='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
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
        elif opt == '--rcutoff':
            corr_cutoff = float(arg)
        elif opt == '--tfcutoff':
            tfcutoff = int(arg)
        elif opt == '--pcutoff':
            pval_cutoff = float(arg)
        elif opt == '--outFile':
            outFile = str(arg)
    
    eExp_df = read_enhancer_exp(eExpCsv)
    
    gExp_df = read_gene_exp(gExpCsv)
    
    merged_exp_df = merge_exp(eExp_df, gExp_df)
    
    eTFBS = read_eTFBS(eTFBS_file, tfcutoff)
    
    gTFBS = read_gTFBS(gTFBS_file, tfcutoff)
    
    enhancers = read_enhancer(eFile)
    
    near_genes = obtain_near_gene(enhancers, gtfFile)
    
    identify_targets(enhancers,
                     near_genes,
                     merged_exp_df,
                     eTFBS,
                     gTFBS,
                     corr_cutoff,
                     pval_cutoff,
                     outFile)