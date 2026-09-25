# Wang-et-al

# Paper Title

Official implementation of:

**"Comparative embryo-mapping reveals neural bias of neuromesodermal progenitors in human stem cell axial elongation models"**

Authors: Yuting Wang; Rafaella Buzatu; Cecile Herbermann; Micha Drukker; Christian Schröter 

## Overview

This repository contains the code and analysis workflows for mapping single-cell RNA-seq data from twelve stem cell-based axial organoid models to stage-matched human and cynomolgus monkey embryos. The analyses characterize progenitor states and developmental trajectories along the human body axis, evaluate how axial organoids recapitulate embryonic development, and investigate the effects of signaling pathway modulation on cell-type composition.


## Installation

Clone the repository:

git clone https://github.com/USERNAME/REPOSITORY.git
cd REPOSITORY

Install dependencies:

pip install -r requirements.txt

## Data

The scRNA-seq datasets are available through the Gene Expression Omnibus (GEO) under accession numbers listed in the paper.
## Requirements

### R

The R-based analyses were performed using R with the following main packages:

- Seurat
- ggplot2
- dplyr
- SingleCellExperiment
- zellkonverter

### Python

The Python-based analyses were performed using Python with the following main packages:

- scanpy
- anndata
- pandas
- numpy
- scikit-learn
- pandas
- scvelo
- loompy
## Analysis workflow

The analysis consists of the following main steps:

1. Preprocessing of scRNA-seq datasets from each oragnoid and human embryo reference
2. Quantification of proportion of each NMP subtype
3. Mapping to human and cynomolgus monkey embryo references
4. Trajectory analysis
5. Multivariate regression analysis
6. Figure generation

## Reproducing the Paper Results

### Figure 1

Run:
`all_mapping.r`

### Figure 2

Run:

`quantify_NMP_subtype.r`

### Figure 3

Run:

`NT_Xue.py`
`NT_Anand.py`
`NT_Sanaki-Matsumiya.py`
`SM_Yamanaka.py`
`SM_Sanaki-Matsumiya.py`
`SM_Miao.py`
`TK_Rito.py`
`TK_Gribaudo.py`
`TK_Makwana.py`
`TK_Sanaki-Matsumiya.py`
`TK_Hamazaki.py`


### Figure 4
Multivariate_regression.py


## License

This project is released under the MIT License.
