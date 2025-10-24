# eModule: A Computational Framework for Enhancer Analysis with Transcriptomic Data

## Overview

eModule is a comprehensive computational framework designed to analyze enhancer expression and identify enhancer-mediated gene regulatory modules using transcriptomic data. This framework addresses the critical need for specialized tools that can effectively leverage RNA-seq data to understand enhancer biology and gene regulation.

## Scientific Background

Enhancers are distal regulatory elements that play crucial roles in gene expression control through multiple mechanisms:
- **Direct transcriptional activation** via transcription factor (TF) binding
- **Chromatin looping** to bring TFs into proximity with target gene promoters
- **Production of enhancer RNAs (eRNAs)** that reflect enhancer activity

eModule specifically focuses on **enhancer-mediated gene modules** - coordinated regulatory units where enhancers facilitate TF-gene interactions through chromatin looping rather than direct TF binding to gene promoters.

## Framework Architecture

eModule consists of three integrated modules that work together to provide comprehensive enhancer analysis:

### 1. Expression Quantification and  Normalization (`calculateExp.py` and `eNormalization.py`)
**Purpose**: Quantifies enhancer expression from RNA-seq bigWig files, and normalizes for cross-sample comparability

**Algorithm**:
- Extracts per-base read coverage across enhancer regions using pyBigWig
- Calculates average coverage values representing eRNA transcriptional activity
- Outputs enhancer loci with corresponding expression values
- Implements TMM (Trimmed Mean of M-values) normalization via edgeR
- Applies CPM (Counts Per Million) transformation
- Filters lowly expressed enhancers (median CPM < 1)
- Produces tissue-specific expression matrices

**Key Features**:
- Processes genome-wide per-base read count data
- Handles multiple samples and tissues efficiently
- Compatible with GTEx bigWig files from recount3 platform
- Integrates R's edgeR package through rpy2 for robust normalization
- Accounts for library size differences and composition bias
- Implements quality filtering to ensure robust downstream analysis

### 2. Differential Expression Analysis (`diffExp.py`)
**Purpose**: Identifies differentially expressed enhancers (e.g., sex-biased)

**Algorithm**:
- Uses generalized linear models (GLM) with Poisson family for count data
- Controls for multiple confounding factors (should be specified in formula file):
  - Demographic variables (age)
  - Technical covariates (RIN, PMI)
  - Population structure (genotype PCs)
  - Gene expression PCs (excluding sex-correlated PCs)
- Applies Benjamini-Hochberg FDR correction
- Filters for biologically significant effects (|log2FC| > 1, FDR < 0.05)

**Key Features**:
- Robust statistical framework appropriate for RNA-seq count data
- Comprehensive covariate adjustment following GTEx consortium methods
- Prevents over-correction by excluding sex-correlated expression PCs

### 3. Regulatory Module Identification (`identifyModule.py`)
**Purpose**: Identifies enhancer-mediated gene regulatory modules

**Algorithm**:
1. **Co-expression Analysis**: Identifies TF-enhancer-gene triplets with significant Spearman correlations
2. **Binding Site Filtering**: Ensures TF binds enhancer but not gene promoter (enhancer-mediated mechanism)
3. **Mediation Analysis**: Uses bootstrap mediation testing to confirm indirect effects
4. **Module Construction**: Groups target genes regulated by the same enhancer-TF pair

**Key Features**:
- Distinguishes enhancer-mediated regulation from co-responsive models
- Implements non-parametric bootstrap mediation analysis
- Uses parallel processing for computational efficiency
- Integrates multiple data types (expression, TFBS, gene annotations)
<p align="center"><img src="https://github.com/LinjieCodes/eModule/blob/main/img/eModule_framework.png" alt="eModule" width="900px" /></p> 

## Installation

### Requirements

**Python Dependencies**:
- Python >= 3.7.0
- numpy >= 1.21.6
- pandas >= 1.3.5
- scipy >= 1.7.3
- statsmodels >= 0.13.5
- pyBigWig >= 0.3.22
- rpy2 >= 3.5.16
- joblib (for parallel processing) >= 1.2.0

**R Dependencies**:
- edgeR >= 3.6.3

**System Requirements**:
- Linux system (tested)
- Sufficient memory for large expression matrices
- Multiple CPU cores recommended for parallel processing

### Installation Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/LinjieCodes/eModule.git
   cd eModule
   ```

2. **Install Python dependencies**:
   ```bash
   pip install numpy pandas scipy statsmodels pyBigWig rpy2 joblib
   ```

3. **Install R dependencies**:
   ```R
   install.packages("BiocManager")
   BiocManager::install("edgeR")
   ```

## Usage

### Data Preparation

**Required Input Files**:

1. **Enhancer Annotations**: BED format file with enhancer coordinates
   - Format: `chromosome\tstart\tend`
   - Example: `chr1\t1000\t2000`

2. **RNA-seq bigWig Files**: Per-base coverage files (e.g., from GTEx/recount3)
   - Can be downloaded from recount3 platform
   - Format: `gtex.base_sums.TISSUE_SAMPLE.1.ALL.bw`

3. **Sample Attributes**: Tab-delimited file with sample metadata
   - Format: `tissue\tsampleID\n`
   - Should include columns for covariates (Sex, Age, RIN, PMI, PCs)

4. **Gene Expression Matrix**: CSV format (genes × samples)
   - Required for regulatory module identification

5. **TFBS Data**: Transcription factor binding site predictions
   - Enhancer TFBS: `enhancer\tTFs\tscore\n`
   - Gene TFBS: `gene\tTFs\tscore\n`
   - Can be generated using JASPAR profiles

6. **Gene Annotations**: GTF file with gene coordinates
   - Required for identifying candidate target genes

### Step-by-Step Analysis

#### Step 1: Quantify Enhancer Expression

```bash
python calculateExp.py \
    --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene \
    --bwfile gtex.base_sums.ADIPOSE_TISSUE_GTEX-1A3MV-2126-SM-718BV.1.ALL.bw \
    --outFile rawExp
```

**Parameters**:
- `--annoFile`: Enhancer annotation file (BED format)
- `--bwfile`: RNA-seq bigWig file
- `--outFile`: Output file for expression values

#### Step 2: Normalize Expression

```bash
python eNormalization.py \
    --annoFile Ensembl_Fantom5_enhancers_nonOverlapGene \
    --sampleFile sampleAttributes \
    --expFolder eRNA_rawExp/ \
    --outFolder eRNA_CPM/
```

**Parameters**:
- `--annoFile`: Enhancer annotation file
- `--sampleFile`: Sample attribute file
- `--expFolder`: Path to folder containing tissue subfolders with expression files. By default, files are named with sample ID. Files of the same tissue are placed in one subfolder.
- `--outFolder`: Output directory for normalized matrices

#### Step 3: Identify Differentially Expressed Enhancers

```bash
python diffExp.py \
    --expFile Spleen_CPM.csv \
    --sampleFile Spleen_sample.csv \
    --tissue Spleen \
	--covarFile  formula.txt \
    --outFile sexBiasedEnhancer
```

**Parameters**:
- `--expFile`: Normalized expression matrix (CSV)
- `--sampleFile`: Sample attributes with covariates
- `--tissue`: Tissue name for output
- `--covarFile`: Formula for GLM fitting
- `--outFile`: Output file for results

#### Step 4: Identify Regulatory Modules

```bash
python identifyModule.py \
    --eExpCsv Spleen_CPM.csv \
    --gExpCsv Spleen_geneExp.csv \
    --eTFBS_file Enhancer_TFBS_demo \
    --gTFBS_file Promoter_TFBS_demo \
    --eFile SpleenSexBiasedEnhancer \
    --gtfFile Homo_sapiens.GRCh38.101.gtf \
    --n_boot 10000 \
    --tfcutoff 400 \
    --rcutoff 0.3 \
    --pcutoff 0.05 \
    --outFile module_out
```

**Parameters**:
- `--eExpCsv`: Enhancer expression matrix
- `--gExpCsv`: Gene expression matrix
- `--eTFBS_file`: Enhancer-TF binding predictions
- `--gTFBS_file`: Gene-TF binding predictions
- `--eFile`: List of enhancers to analyze
- `--gtfFile`: Gene annotation file
- `--n_boot`: Number of bootstrap resamples (default: 10000)
- `--tfcutoff`: TF binding score threshold (default: 400)
- `--rcutoff`: Correlation coefficient threshold (default: 0.3)
- `--pcutoff`: P-value threshold (default: 0.05)

## Output Formats

### Expression Quantification Output
```
locus	average_coverage
chr1:1000-2000	15.67
chr1:5000-6000	8.23
...
```

### Differential Expression Output
```
Tissue	Enhancer	Male_exp	Female_exp	log2(FC)	coef	Pval	FDR
Spleen	chr1:1000-2000	4.567	2.123	2.444	0.567	1.2e-05	0.0012
...
```

### Regulatory Module Output
```
Enhancer	Regulating TF	Target genes
chr1:1000-2000	FOXA1	GENE1, GENE2, GENE3
chr1:5000-6000	GATA3	GENE4, GENE5
...
```

## Biological Interpretation

### Enhancer-Mediated vs. Co-responsive Models

eModule specifically identifies **enhancer-mediated regulation** by requiring:
1. TF binding to enhancer (high confidence TFBS)
2. TF NOT binding to gene promoter (no TFBS)
3. Significant co-expression among TF-enhancer-gene triplet
4. Significant mediation effect (indirect pathway)

This approach distinguishes from alternative models:
- **Co-responsive model 1**: All three components respond to upstream regulator
- **Co-responsive model 2**: TF binds both enhancer and gene promoter

### Validation Strategies

eModule results can be validated through:
1. **Chromatin interaction data** (Hi-C, ChIA-PET)
2. **eQTL analysis** (GTEx cis-eQTLs)
3. **Functional enrichment** (GO terms, pathways)
4. **Expression divergence** between conditions

## Example Applications

### Sex-Biased Enhancer Analysis
- Identify enhancers with sex-specific activity
- Understand mechanisms of sex-biased gene expression
- Link enhancer variation to phenotypic differences

### Disease-Associated Enhancer Studies
- Compare enhancer activity between disease and control samples
- Identify dysregulated enhancer-gene modules
- Understand enhancer contribution to disease mechanisms

### Developmental Enhancer Dynamics
- Track enhancer activity changes during development
- Identify stage-specific regulatory modules
- Understand enhancer roles in cell fate decisions

## Performance Considerations

### Computational Requirements
- **Memory**: Proportional to number of enhancers × samples
- **Time**: Module identification scales with number of TF-enhancer-gene triplets
- **Parallelization**: Bootstrap mediation uses joblib for parallel processing

### Optimization Tips
1. **Filter enhancers** before analysis to reduce computational burden
2. **Use appropriate correlation thresholds** to limit triplet combinations
3. **Parallel processing** automatically utilizes available CPU cores
4. **Memory-efficient** data structures minimize RAM usage

## Troubleshooting

### Common Issues

1. **Import errors with rpy2**:
   - Ensure R is properly installed and in PATH
   - Check R_HOME environment variable
   - Verify edgeR package installation

2. **Memory errors with large datasets**:
   - Process tissues separately
   - Filter lowly expressed enhancers earlier
   - Use smaller bootstrap sample sizes

3. **Slow performance**:
   - Increase correlation threshold to limit triplet combinations
   - Use fewer bootstrap resamples for initial testing
   - Process chromosomes or regions separately

### Data Quality Checks

1. **Verify bigWig file integrity**:
   ```python
   import pyBigWig
   bw = pyBigWig.open("file.bw")
   print(bw.header())
   ```

2. **Check expression matrix dimensions**:
   ```python
   import pandas as pd
   df = pd.read_csv("expression.csv", index_col=0)
   print(df.shape)
   print(df.head())
   ```

3. **Validate TFBS file format**:
   ```bash
   head -n 5 Enhancer_TFBS.txt
   ```

## Citation

If you use eModule in your research, please cite:

Jie Lin et al. *xxx*, xxx, xxx.

## Contributing

We welcome contributions to eModule:
- Reporting bugs
- Suggesting features
- Submitting pull requests
- Improving documentation

## License

This project is licensed under the Apache 2.0 License - see LICENSE file for details.

## Contact

For questions, bug reports, or collaboration inquiries:
- Email: JL.linjie@outlook.com
- GitHub Issues: https://github.com/LinjieCodes/eModule/issues

## Acknowledgments

- **GTEx Consortium** for providing transcriptomic data
- **recount3 project** for processed bigWig files
- **FANTOM and ENCODE** projects for enhancer annotations
- **JASPAR** database for TF binding profiles
- **edgeR developers** for robust normalization methods
- **pyBigWig developers** for a python extension for quick access to bigBed files