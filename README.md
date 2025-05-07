# ENSAS - Expression Neighborhood Sequence Association Study
Identifying associations of de novo noncoding variants with phenotypes through integration of gene expression, sequence and sex information
https://www.biorxiv.org/content/10.1101/2024.03.20.585624v1

## Requirements
```
numpy==1.24.3
pandas==2.0.3
scikit_learn==1.3.0
scipy==1.12.0
```

## Basic usage
```
usage: model.py [-h] [-v [VARIANTS]] [-g [GENES]] [-o [OUTPUT]]
                [-n [NEIGHBORS]] [-k [KMER]] [-c [GENE_COL]]

options:
  -h, --help            show this help message and exit
  -v [VARIANTS], --variants [VARIANTS]
                        Variants file
  -g [GENES], --genes [GENES]
                        Gene matrix file
  -o [OUTPUT], --output [OUTPUT]
                        Output path
  -n [NEIGHBORS], --neighbors [NEIGHBORS]
                        Size of each neighborhood, default 1000
  -k [KMER], --kmer [KMER]
                        Length of k-mers, default 6
  -c [GENE_COL], --gene_col [GENE_COL]
                        Column for gene assignment, default: GeneNAMEOuterTSS
```

The input variants file must have columns "Proband" (binary proband-sibling phenotype labels), "sequence" (DNA sequence), and a column specifying gene assignment for each variant (default: "GeneNAMEOuterTSS"). Sample files can be found in the ```data``` folder. 
The gene matrix file must be tab-delimited, with a header and an index column. A sample gene matrix file can be found in the ```data``` folder.

## Output format
The columns in the output file are: 
- Gene name
- The local GC content Mann-Whitney U effect size, computed from dividing the U-statistic by the total number of comparisons
- The local GC content Mann-Whitney U p-value
- The K-mer Naive Bayes model Mann-Whitney U effect size on the testing fold
- The K-mer Naive Bayes model Mann-Whitney U p-value on the testing fold
- The local GC content Mann-Whitney U effect size on the testing fold
- The local GC content Mann-Whitney U p-value on the testing fold
All p-values reported are uncorrected. 
