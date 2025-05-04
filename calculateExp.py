import os
import sys, getopt
#import multiprocessing
import pyBigWig


def read_enhancer(annoFile):
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
    bw = pyBigWig.open(bwfile)
    counts = {}
    for chrom, start, end, locus in enhancers:
        count = bw.stats(chrom, start, end, exact=True)[0] #the average value over a range
        counts[locus] = str(round(count, 2))
    bw.close()
	
    with open(outFile, 'w') as f_re:
        for t in enhancers:
            locus = t[-1]
            count = counts[locus]
            f_re.write('\t'.join([locus, count])+'\n')


def usage():
    print("""Parameters:
        --annoFile       enhancer_annotation_file. 
        --bwfile         bigwig file. This file can be downloaded from the Recount3 platform.
        --outFile        the file to write result.
        """)
    print()
    print('Example: python calculateExp.py  --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene --bwfile gtex.base_sums.ADIPOSE_TISSUE_GTEX-1A3MV-2126-SM-718BV.1.ALL.bw --outFile rawExp')
    
       
if __name__ == '__main__':
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['annoFile=', 'bwfile=', 'outFile='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
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
            
    enhancers = read_enhancer(annoFile)
    
    extract_exp(bwfile, enhancers, outFile)