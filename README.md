# aging-gtex-xci

#### Repository for an analysis of age-related escape from XCI

This repository contains scripts used in the analysis of XCI in three females from GTEx with nmXCI.

I ran these analyses using [R](https://cran.r-project.org/) (v4.3).

# Inputs

The following data files are expected:

* GTEx nmXCI ASE data from Gylemo et al: gylemo.csv
* TADS: A-172_GSE147123_tad.bed
  
# Pipeline

### all scripts required for analysis and figures

```
# Create sex-specific splici transcriptomes and references
scripts/Rscripts.R
```
