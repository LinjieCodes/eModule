# A computational framework (enhancer network, eNet) for enhancer analysis with transcritomic data

The eNet program comprises four distinct scripts, each serving a specific purpose. These purposes are as follows: firstly, to quantify enhancer expression levels; secondly, to conduct expression normalization; thirdly, to identify differentially expressed enhancers; and lastly, to infer enhancer-mediated gene regulation networks. Collectively, these scripts enable the program to analyze enhancer-related gene regulation with transcritomic data.

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **rpy2**: >=3.5.16

6. **edgeR**: >=3.6.3

7. **statsmodels**: >=0.13.5

8. **pyBigWig**: >=0.3.22

9. **OS**: the eNet code has been tested on Linux system.


#1. Calculate expression of enhancer

## Data
1. **Ensembl_Fantom5_enhancers_nonOverlapGene**  --  the annotation file of all enhancers.

2. **gtex.base_sums.ADIPOSE_TISSUE_GTEX-1A3MV-2126-SM-718BV.1.ALL.bw**  --  the bigwig file downloaded from the Recount3 platform.


## Usage
Here is a command line to run the calculateExp.py script:

```
'Example: python calculateExp.py  --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene --bwfile gtex.base_sums.ADIPOSE_TISSUE_GTEX-1A3MV-2126-SM-718BV.1.ALL.bw --outFile rawExp'
```

## Help information
Here is a brief explanation of the command line arguments:

```
Parameters       Functions
--annoFile       enhancer_annotation_file. 
--bwfile         bigwig file. This file can be downloaded from the Recount3 platform. 
--outFile        the file to write result.
```


#2. Normalization of enhancer expression

## Data
1. **Ensembl_Fantom5_enhancers_nonOverlapGene**  --  the annotation file of all enhancers.

2. **sampleAttributes**  --  the sample attribute file.

3. **eRNA_perbase_average**  --  the folder of enhancers' raw expression.

## Usage
Here is a command line to run the eNormalization.py script:

```
'Example: python eNormalization.py  --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene --sampleFile sampleAttributes --expFolder eRNA_perbase_average/ --outFolder eRNA_RPM/'
```

## Help information
Here is a brief explanation of the command line arguments:

```
Parameters       Functions
--annoFile       enhancer_annotation_file. 
--sampleFile     sample_attribute_file. 
--expFolder      Path of enhancers' raw expression.
--outFolder      Output path to write enhancers' normalized expression.
```


#Identification of differentially expressed enhancers

## Data
1. **Spleen_RPM.csv**  --  the enhancer expression matrix.

2. **Spleen_sample.csv**  --  sample_attribute_file of a tissue. This file can be downloaded from the Recount3 platform.

## Usage
A command line to run the diffExp.py script:

```
'Example: python diffExp.py  --expFile Spleen_RPM.csv --sampleFile Spleen_sample.csv --tissue Spleen --outFile sexBiasedEnhancer'
```

## Help information
A brief explanation of the command line arguments:

```
Parameters      Functions
--expFile       enhancer RPM matrix. 
--sampleFile    sample_attribute_file of the tissue. 
--tissue        the tissue label.
--outFile       the file to write result.
```

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com.
