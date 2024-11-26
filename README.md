# ctDRTF
A network-based polygenic enrichment method for identifying cell type-specific transcription factor-regulatory programs relevant to complex diseases.

# ================================
@ Latest version: v2.0.0 (2024-11-24: scMORE)
@ Previous Versions:
    - v1.2.0 (2024-04-15: ctDRTF)
    - v1.1.1 (2024-03-15: ctDRTF)
    - v1.1.0 (2024-01-01: ctDRTF)
    - v1.0.0 (2023-11-03: ctDRTF)
# -------------------------------

![logo](https://github.com/mayunlong89/ctDRTF/blob/main/example/figure/Picture1.png)

ctDRTF assumes that if GWAS-identified disease-specific genes are concordantly activated in a cell type-specific TF-regulon, then the TF is more likely to have a pivotal role in disease via the given cell type. Thus, the computational framework ctDRTF (cell type-specific Disease-Relevant Transcription Factor) is designed for performing regulatory network-based inference of the associations between TF-related regulons and disease-specific gene sets in a context-specific manner. The input of ctDRTF includes multimodal matrix of both scRNA-seq and scATAC-seq data, and GWAS summary data for a quantitative trait or disease (case-control study). 

![Workflow](https://github.com/mayunlong89/ctDRTF/blob/main/example/figure/Figure_1_v3.png)


# Installing ctDRTF
We recommend installing ctDRTF via github using devtools:

```r

library(devtools)
install_github("mayunlong89/ctDRTF")

```
See the DESCRIPTION file for a complete list of R dependencies. If the R dependencies are already installed, installation should finish in a few minutes.

## How to run ctDRTF
```r
#@' main function
#@' single_cell: the input single-cell data.
#@' n_genes: the minimum number of genes in a given regulon 
#@' MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
#@' MAGMA_GWAS_data header: SYMBOL, logP, ZSTAT
#@' Gene_num: The number of disease-specific genes, default set to 500
#@' MC_num: Set 100 times of Monte Carlo simulation
#@' theta range from 0.1 ~ 1, default set to 0.5
#@' mi (power, or beta): expand the specificity difference between cell types,default set to 1
#@' mode: mo=1 indicates use magma-based z-scores as weights; mo=0 indicates no weights

ctdrtf <- function(single_cell = single_cell,
                   MAGMA_GWAS_data = MAGMA_GWAS_data,
                   n_genes= 10,
                   Gene_num = 500,
                   MC_num = 1000,
                   theta=0.5,
                   mi=1,
                   mo=1)

```


| Function                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `ctdrtf()`                         | The main function to fit the model for identification of cell type-specific regulons relevant to complex diseases           |
| `GRN_func()`                       | Construct global TF-gene regulotory network depended on the Pando framework                                                 |
| `COSR_pre_func()`                  | Calculate the cell type-level specificity score of each gene or TF using the cosine similarity algorithm                    |
| `COSR_func_weight()`               | Identify cell type-specific regulons relevant to disease using polygenic enrichment method                                  |
| `MC_JSI_score_func_weight()`       | Calculate the empirical P value for each regulon-disease association using Monote Carlo permutation algorithm               |
| `max_min_scale()`                  | Scale the specificity score between 0 and 1 across all cells                                                                |




### The ctDRTF pipeline             
```r
# Step 1
##@ 1) Constructing global TF-gene regulatory network
##single_cell: the input single-cell data--Seurat object.
##n_genes: the minimum number of genes in a given regulon`

GRN_func(single_cell = single_cell,n_genes=n_genes)

# Step 2
##@ 2) Calculating specificity scores
##single_cell is the input single-cell data.
##num_genes <- length(rownames(single_cell))

COSR_pre_func(single_cell=single_cell)

# Step 3
@' 3) Identifying cell type-specific regulons relevant to disease
#@' magma result processing
#@' magma_results <- magma_results %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
#@' MAGMA_GWAS_data <- magma_results[,c(10,11,8)]
#@' MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
#@' MAGMA_GWAS_data header: SYMBOL, logP, ZSTAT
#@' data_regulons1: TF-gene pairs matrix
#@' tf_left: all tf names
#@' MC_num: Set 1000 times of running MC_JSI_score_func()
#@' data_s1: matrix of genes and TFs specificity scores across cell types
#@' theta range from 0.1 ~ 1, default set to 0.5
#@' mi(power): expand the specificity difference between cell types, default set to 1
#@' mode: mo=1 indicates use magma-based z-scores as weights; mo=0 indicates no weights

COSR_func_weight <- function(tf_left=tf_left,
                      data_s1=data_s1,
                      data_regulons1=data_regulons1,
                      MAGMA_GWAS_data = MAGMA_GWAS_data,
                      MC_num = 1000,
                      Gene_num=500,
                      theta=0.5,
                      mi=1,
                      mo=1)


# Step 4
#@ 4) Calculating the empirical P value for each regulon-disease association
#@' MC Function
#@' data_s1_sub: The input of matrix containing the specificity score of genes in each cell type
#@' tf_left: All the high-qualified TFs that passed the quality control
#@' tf_left_1 <- tf_left[which(tf_left!=tf_left[j])] #removing the targeted TF as controls
#@' len_of_regulon = length(M1_regulon), the number of genes and TFs in a given regulon
#@' all_genes: All genes from phenotype-associations based on MAGMA,ie.,length(MAGMA_GWAS_data$SYMBOL)
#@' MAGMA_GWAS_data header: SYMBOL, logP, ZSTAT
#@' Gene_num: The number of disease-specific genes, default set to 500
#@' theta range from 0.1 ~ 1, default set to 0.5
#@' mi (power): expand the specificity difference between cell types,default set to 1
#@' mode: mo=1 indicates use magma-based z-scores as weights; mo=0 indicates no weights
#@' Random-based specificity*JSI score for each regulon

MC_JSI_score_func_weight<- function(data_s1_sub = data_s1_sub,
                                    tf_left_1 = tf_left_1,
                                    len_of_regulon = len_of_regulon,
                                    all_genes = all_genes, 
                                    Gene_num = 500,
                                    theta=0.5,
                                    mi=1,
                                    mo=1)

```

## Assigning cell types to single-cell data

```r

Idents(single_cell) <- single_cell$cell_type

```

### Generate MAGMA-based gene set

```shell
1) MAGMA codes for generating disease-relevant genes

#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/MAGMA
export DATA=/share/pub/mayl/MAGMA_test
export OUTPUT=/share/pub/mayl/MAGMA_test

#MAGMA annotation:

$MAGMA_DIR/magma \
    --snp-loc  $DATA/GWAS_UKBiobank_summary_final.hg19.location  \
    --annotate window=20,20 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation  

#gene-based association analysi:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/GWAS_UKBiobank_summary_final.results_Pval \
    N=13239 \
    --gene-annot   $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation.genes.annot  \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P


2) Processing MAGMA-results: 'magma.genes.out'
#MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
#header of MAGMA_GWAS_data: SYMBOL, logP, ZSTAT

magma_results <- read.table("magma.genes.out",header = TRUE)
magma_results <- magma_results %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
MAGMA_GWAS_data <- magma_results[,c(10,11,8)]


```
#For more detailed codes on MAGMA tool, please refer to [here](https://cloufield.github.io/GWASTutorial/09_Gene_based_analysis/)

### Example format

```
1) Single-cell data

The input format of single-cell data: Seurat-generated S4 object.

ctDRTF is fully compatiable with Seurat, a widely-used single-cell analysis tool.


2) MAGMA-result data (i.e., MAGMA_GWAS_data)

  SYMBOL      logP  ZSTAT
1     PON1 10.899319 6.6721
2     PON2  8.570943 5.8352
3     PON3  4.386613 3.9381
4 CDKN2AIP  4.130850 3.7944
5     CASR  4.083372 3.7672
6     ASB4  3.919446 3.6719

````

# Citations
1. Ma et al., Polygenic network enrichment identifies cellular context-specific regulons relevant to diseases by integration of single-cell multiomic data, `Genome Biology`,(under review), 2024
2. Ma et al., Sytematic dissection of pleiotropic loci and critical regulons in exhibitory neurons and microglia relevant to neuropsychiatric and ocular diseases, [Research Square](https://www.researchsquare.com/article/rs-4514542/v1), 2024.


# scHOB Database
Human organoids are advanced three-dimensional structures that accurately recapitulate key characteristics of human organ development and functions. Unlike two-dimensional cultures lacking critical cell-cell communications, organoids provides a powerful model for recovering complex cellular dynamics involved in developmental and homeostatic processes. Organoids also allow genetic and pharmacological manipulation in a more physiologically relevant context compared to animal models. Although single-cell sequencing advancements have accelerated their biological and therapeutic use, there has been no systematic platform for unified processing and analysis of organoid-based single-cell multiomics data. 

We thus established scHOB (single-cell Human Organoid Bank), a multi-omic single-cell database, consisting of both scRNA-seq and scATAC-seq data on 10 types of widely-adopted human organoids (i.e., brain, lung, heart, eye, liver & bile duct, pancreas, intestine, kidney, and skin) spanning more than 1.5 million cells with 67 main cell types in 385 samples across 83 distinct protocols. see [Github code](https://github.com/mayunlong89/scHOB/tree/main); see [scHOB Website](https://schob.su-lab.org/).
The single-cell multiome data in scHOB have been used by ctDRTF, see Ma et al. 2024.

# Application example of scHOB database:
Ma et al., Integration of human organoids single-cell transcriptomic profiles and human genetics repurposes critical cell type-specific drug targets for severe COVID-19. [Cell Proliferation](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558),2024, and see related [Github codes](https://github.com/mayunlong89/scHuman_organoids_COVID19).


# Other references:
 
1. [scPagwas](https://www.cell.com/cell-genomics/pdf/S2666-979X(23)00180-5.pdf)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8137370.svg)](https://doi.org/10.5281/zenodo.8137370)

2. Development of novel polygenic regression method scPagwas for integrating scRNA-seq data with GWAS on complex diseases. see [Ma et al. Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5), and see related [Github codes](https://github.com/mayunlong89/scPagwas_main)









