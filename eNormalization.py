import os
import sys, getopt
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri
r.source('tmm.r')

pandas2ri.activate()


#read enhancer annotation
def read_enhancer(annoFile):
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


#read tissues and samples
def read_sample(sampleFile):
    tissue_samples = {}
    all_samples = set()
    with open(sampleFile) as f:
        f.readline()
        for line in f:
            cols = line.split('\t')
            tissue = cols[0]
            sample = cols[1]
            if tissue not in tissue_samples:
                tissue_samples[tissue] = []
            tissue_samples[tissue].append(sample)
            all_samples.add(sample)
    return tissue_samples, all_samples
    

#read eRNA expression
def read_enhancerExp(enhancer_index, 
                     lengths, 
                     all_samples, 
                     expFolder):
    sample_counts = {}
    for tissue in os.listdir(expFolder):
        print(tissue)
        for sampleFile in os.listdir(expFolder+tissue):
            sampleID = sampleFile[sampleFile.find('GTEX-'):(sampleFile.find('.ALL.bw')-2)]
            if sampleID in all_samples:
                sample_counts[sampleID] = []
                this_sample_counts = {}
                with open(expFolder+tissue+'/'+sampleFile) as f:
                    for line in f:
                        cols = line.split('\t')
                        if cols[0] in lengths:
                            this_sample_counts[cols[0]] = round(float(cols[1]) * lengths[cols[0]])
                for enhancer in enhancer_index:
                    sample_counts[sampleID].append(this_sample_counts[enhancer])
    return sample_counts


#enhancer expression normalization
def normalize_exp(enhancer_index,
                  tissue_samples,
                  sample_counts,
                  outFolder):    
    #convert to dataframe with enhancer as index, and sample as column
    sampleList = list(sample_counts.keys())
    expValues = []
    for sample in sampleList:
        expValues.append(sample_counts[sample])
    counts_df = pd.DataFrame(expValues, columns=enhancer_index, index=sampleList)
    counts_df = counts_df.T

    #convert to R dataframe
    counts_df_r = pandas2ri.py2rpy(counts_df)

    #TMM normalization
    outFile = 'TMM_read_counts.csv'
    r.tmm_nor(counts_df_r, outFile)

    outFile = 'TMM_read_counts.csv'
    #read TMM-normalized read counts
    counts_tmm = pd.read_csv(outFile, index_col=0)

    #normalized by reads per million (RPM) method
    constant = 1000000
    counts_rpm = counts_tmm.apply(lambda x: round((x / x.sum()) * constant))
    counts_rpm = counts_rpm.astype(int)

    #export csv file for each tissue
    for tissue in tissue_samples:
        samples = tissue_samples[tissue]
        
        counts_rpm_tissue = counts_rpm[samples]
        
        medians = counts_rpm_tissue.median(axis=1)
        
        columns_to_drop = medians[medians < 1].index

        counts_rpm_tissue = counts_rpm_tissue.drop(index=columns_to_drop)
        
        counts_rpm_tissue.to_csv(outFolder + '%s.csv' % tissue, index=True)    


def usage():
    print("""Parameters:
        --annoFile       enhancer_annotation_file. 
        --sampleFile     sample_attribute_file. 
        --expFolder      Path of enhancers' raw expression.
        --outFolder      Output path to write enhancers' normalized expression.
        """)
    print()
    print('Example: python eNormalization.py  --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene --sampleFile sampleAttributes --expFolder eRNA_perbase_average/ --outFolder eRNA_RPM/')

    
if __name__ == '__main__':
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['annoFile=', 'sampleFile=', 'expFolder=', 'outFolder='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
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
     
    if not os.path.exists(outFolder):
        os.mkdir(outFolder)
    
    enhancer_index, lengths = read_enhancer(annoFile)

    tissue_samples, all_samples = read_sample(sampleFile)
    
    sample_counts = read_enhancerExp(enhancer_index, lengths, all_samples, expFolder)

    normalize_exp(enhancer_index,
                  tissue_samples,
                  sample_counts,
                  outFolder)  