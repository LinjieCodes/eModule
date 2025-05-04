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


#identify differentally expressed enhancers
def diff_exp(expFile,
             sampleFile,
             tissue,
             outFile):
    results = []
    p_values = []

    # read expression data
    expression_data = pd.read_csv(expFile, index_col=0)
    
    # Samples as index, enhancers as column
    expression_data = expression_data.T
    
    # read attribute data
    individual_attributes = pd.read_csv(sampleFile, index_col=0)
    
    # Samples as index, attribute as column
    individual_attributes = individual_attributes.T
    
    #identify gene expression PCs that significantly correlate with sex by Point-Biserial Correlation test
    correlatedPCs = set()
    expPCs = [PC for PC in individual_attributes.columns if 'InferredCov' in PC]
    for expPC in expPCs:
        corr, p_value = pointbiserialr(individual_attributes['Sex'], individual_attributes[expPC])
        if p_value < 0.05:
            correlatedPCs.add(expPC)
            
    #exclude sex-correlated gene expression PCs
    for expPC in correlatedPCs:
        individual_attributes = individual_attributes.drop(expPC, axis = 1)
    
    # construct a formula for GLM
    formula = 'expression ~ ' + ' + '.join(individual_attributes.columns)
    
    # Iterate over each enhancer
    for enhancer in expression_data.columns:
        enhancer_expression_data = expression_data[[enhancer]]
        
        #enhancer_expression_data.rename(columns = {enhancer: 'expression'}, inplace=True)
        enhancer_expression_data.columns  = ['expression']
        
        # median expression for male samples
        male_median = round(float(np.median(np.log2(enhancer_expression_data[individual_attributes['Sex']==1]+0.01))), 3)
        
        # median expression for female samples
        female_median = round(float(np.median(np.log2(enhancer_expression_data[individual_attributes['Sex']==2]+0.01))), 3)
        
        # male VS female fold change
        fold_change = round(male_median - female_median, 3)
    
        combined_data = enhancer_expression_data.join(individual_attributes)
    
        # Fit a generalized linear model
        model = glm(formula, combined_data, family=Poisson()).fit()

        # identify sex-associated enhancers
        # extract p values
        model_summary = model.summary2().tables[1]
        pval = model_summary['P>|z|']['Sex']
        coef = round(model_summary['Coef.']['Sex'], 3)
        
        p_values.append(pval)
        results.append([tissue, enhancer, male_median, female_median, fold_change, coef, format(pval, '.3e')])

    #Benjamini-Hochberg p-value correction
    p_values = np.array(p_values)
    r, fdrs, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(p_values, alpha=0.05, method='fdr_bh')
    for i in range(len(fdrs)):
        results[i].append(fdrs[i])

    #write result                         
    with open(outFile, 'w') as f:
        f.write('\t'.join(['Tissue',
                           'Enhancer',
                           'Male_exp',
                           'Female_exp',
                           'log2(FC)',
                           'coef',
                           'Pval',
                           'FDR'])+'\n')
        for tissue, enhancer, male_median, female_median, fold_change, coef, pval, fdr in results:
            if fdr < 0.05 and abs(fold_change)>1 and coef*fold_change<0:
                f.write('\t'.join([tissue,
                                   enhancer,
                                   str(male_median),
                                   str(female_median),
                                   str(fold_change),
                                   str(coef),
                                   str(pval),
                                   str(format(fdr, '.3e'))])+'\n')    


def usage():
    print("""Parameters:
        --expFile        enhancer RPM matrix. 
        --sampleFile     sample_attribute_file of the tissue. The file can be downloaded from Recount3 platform.  
        --tissue         the tissue label.
        --outFile        the file to write result.
        """)
    print()
    print('Example: python diffExp.py  --expFile Spleen_RPM.csv --sampleFile Spleen_sample.csv --tissue Spleen --outFile sexBiasedEnhancer')
    

if __name__ == '__main__':
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['expFile=', 'sampleFile=', 'tissue=', 'outFile='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
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
    
    diff_exp(expFile,
             sampleFile,
             tissue,
             outFile)
