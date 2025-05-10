# A computational framework (eModule) for enhancer analysis with transcritomic data

The eNet program comprises four scripts, each serving a specific purpose. These purposes are as follows: firstly, to quantify enhancer expression levels; secondly, to conduct enhancer expression normalization; thirdly, to identify differentially expressed enhancers; and lastly, to infer enhancer-mediated gene regulation modules. Collectively, these scripts enable the program to analyze enhancer-related gene regulation with transcritomic data.

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


# 1. Calculate expression of enhancer

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
--annoFile       to specify enhancer_annotation_file. 
--bwfile         bigwig file. This file can be downloaded from the Recount3 platform. 
--outFile        the file to write result.
```


# 2. Normalization of enhancer expression

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
--annoFile       to specify enhancer annotation file. 
--sampleFile     to specify sample attribute file. 
--expFolder      the folder of enhancers' raw expression.
--outFolder      to specify the utput folder to write enhancers' normalized expression.
```


# 3.Identification of differentially expressed enhancers

## Data
1. **Spleen_RPM.csv**  --  enhancer expression matrix.

2. **Spleen_sample.csv**  --  sample attribute. This file can be downloaded from the Recount3 platform.

## Usage
A command line to run the diffExp.py script:

```
'Example: python diffExp.py  --expFile Spleen_RPM.csv --sampleFile Spleen_sample.csv --tissue Spleen --outFile sexBiasedEnhancer'
```

## Help information
A brief explanation of the command line arguments:

```
Parameters      Functions
--expFile       to specify enhancer RPM matrix. 
--sampleFile    to specify sample_attribute_file. 
--tissue        the tissue label.
--outFile       the file to write result.
```


# 4.Inference of enhancer-mediated gene regulation modules

## Data
1. **Spleen_RPM.csv**  --  enhancer expression matrix.

2. **Spleen_geneExp.csv**  --  gene expression matrix.

3. **Enhancer_TFBS_demo**  --  Enhancers' binding TFs. This is a demo showing the file format. The complete file is too large to upload, which can be obtained through the JASPAR website.

4. **Promoter_TFBS_demo**  --  Genes' binding TFs. This is a demo to showing the file format. The complete file is too large to upload, which can be obtained through the JASPAR website.

5. **SpleenSexBiasedEnhancer**  --  List of enhancers you are interested.

6.  **Homo_sapiens.GRCh38.101.gtf**  -- The gtf file can be downloaded from the Ensembl website. This file is too large to upload.

## Usage
A command line to run the identifyModule.py script:

```
'Example: python identifyModule.py  --eExpCsv Spleen_RPM.csv  --gExpCsv Spleen_geneExp.csv --eTFBS_file Enhancer_TFBS_demo --gTFBS_file Promoter_TFBS_demo --eFile SpleenSexBiasedEnhancer --gtfFile Homo_sapiens.GRCh38.101.gtf --tfcutoff 400 --rcutoff 0.3 --pcutoff 0.05 --outFile module_out'
```

## Help information
A brief explanation of the command line arguments:

```
Parameters      Functions
--eExpCsv       to specify enhancer expression matrix. 
--gExpCsv       to specify gene expression matrix. 
--eTFBS_file    to specify the file of enhancer-TF binding pairs.
--gTFBS_file    to specify the file of gene-TF binding pairs.
--eFile         list of enhancers to be analyzed.
--gtfFile       the gtf file.
--tfcutoff      the score cutoff to determine potential enhancer-TF/gene-TF bindings.
--rcutoff       correlation coefficient cutoff for Spearman correlation test.
--pcutoff       p-value cutoff for Spearman correlation test.
--outFile       to specify the output file.
```

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com.
