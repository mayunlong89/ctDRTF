# ctDRTF
A computational method for identifying cell type-specific transcription factor relevant to complex disease.

![logo](https://github.com/mayunlong89/ctDRTF/blob/main/example/figure/Picture2.png)

ctDRTF assumes that if GWAS-identified disease-specific genes are concordantly activated in a cell type-specific TF-regulon, then the TF is more likely to have a pivotal role in disease via the given cell type. Thus, the computational framework ctDRTF (cell type-specific Disease-Relevant Transcription Factor) is designed for performing regulatory network-based inference of the associations between TF-related regulons and disease-specific gene sets in a context-specific manner. The input of scDRTF includes multimodal matrix of both scRNA-seq and scATAC-seq data, and GWAS summary data for a quantitative trait or disease (case-control study). 

![Workflow](https://github.com/mayunlong89/ctDRTF/blob/main/example/figure/ctDRTF_workflow_only.png)


# Installing ctDRTF
We recommend installing ctDRTF via github using devtools:

```r

library(devtools)
install_github("mayunlong89/ctDRTF")

```
See the DESCRIPTION file for a complete list of R dependencies. If the R dependencies are already installed, installation should finish in a few minutes.

## How to run ctDRTF
#main function
#single_cell: the input single-cell data.
#n_genes: the minimum number of genes in a given regulon.
#MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P).
#Gene_num: The number of disease-specific genes, default set to 500.
#MC_num: Set 100 times of running MC_JSI_score_func().
#theta range from 0 ~ 1, default set to 0.5.
#mode: default "weight", alternatively, "none"; This parameter is used the z-score of each gene from magma as weight.

```r

ctdrtf <- function(single_cell = single_cell,
                           MAGMA_GWAS_data = MAGMA_GWAS_data,
                           n_genes= 10,
                           Gene_num = 500,
                           MC_num = 1000,
                           theta=0.5,
                           mode="weight")

```


| Function                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `sccomp_estimate`                  | Fit the model onto the data, and estimate the coefficients                                                                  |
| `sccomp_remove_outliers`           | Identify outliers probabilistically based on the model fit, and exclude them from the estimation                            |
| `sccomp_test`                      | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate`                 | Simulate data from the model, or part of the model                                                                          |
| `sccomp_predict`                   | Predicts proportions, based on the mode, or part of the model                                                               |
| `sccomp_remove_unwanted_variation` | Removes the variability for unwanted factors                                                                                |
| `plot`                             | Plors summary plots to asses significance                                                                                   |



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
##@ 3) Identifying cell type-specific regulons relevant to disease
##MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
##data_regulons1: TF-gene pairs matrix
##tf_left: all tf names
##MC_num: Set 1000 times of running MC_JSI_score_func()
##data_s1: matrix of genes and TFs specificity scores across cell types
##theta range from 0 ~ 1, default set to 0.5

#mode 1
COSR_func_weight <- function(tf_left=tf_left,
                      data_s1=data_s1,
                      data_regulons1=data_regulons1,
                      MAGMA_GWAS_data = MAGMA_GWAS_data,
                      MC_num = 1000,
                      Gene_num=500,
                      theta=0.5)

#mode 2
COSR_func(tf_left=tf_left,
          data_s1=data_s1,
          data_regulons1=data_regulons1,
          MAGMA_GWAS_data = MAGMA_GWAS_data,
          MC_num = 1000,
          Gene_num=500,
          theta=0.5)

# Step 4
#@ 4) Calculating the P value for each regulon-disease link
#MC Function
#data_s1_sub: The input of matrix containing the specificity score of genes in each cell type
#tf_left: All the high-qualified TFs that passed the quality control
#tf_left_1 <- tf_left[which(tf_left!=tf_left[j])] #removing the targeted TF as controls
#len_of_regulon = length(M1_regulon), the number of genes and TFs in a given regulon
#all_genes: All genes from phenotype-associations based on MAGMA,ie.,length(MAGMA_GWAS_data$SYMBOL)
#Gene_num: The number of disease-specific genes, default set to 500

MC_JSI_score_func(data_s1_sub = data_s1_sub,
                             tf_left_1 = tf_left_1,
                             len_of_regulon = len_of_regulon,
                             all_genes = all_genes, 
                             Gene_num = 500)

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

### Citation
Ma et al., Unveiling cell type-specific transcriptional regulomes from single-cell transcriptomic and epigenomic atlas of human organoids, medRvix, 2024












