
# How to run example using ctDRTF

```
#Single-cell multi-omic data on pbmc_10x (Illumina Official Website)
single_cell <- readRDS("pbmc_10x_example.rds")

#MAGMA-based gene-based analysis of GWAS on monocyte count (UKBiobank ID:ieu-b-31)
magma_example <- readRDS("magma_example.rds")

#running ctDRTF
data_example <- ctdf_main_func(single_cell = single_cell,
                               MAGMA_GWAS_data = magma_example,
                               n_genes = 10,
                               Gene_num = 500,
                               MC_num=100,
                               theta=0.5)

```


## Assigning cell types to single-cell data

```
Idents(single_cell) <- single_cell$cell_type

```


### Generate MAGMA-based gene set

```
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
magma_results <- read.table("magma.genes.out",header = TRUE)
magma_results <- magma_results %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
MAGMA_GWAS_data <- magma_results[,c(10,11)]


```
#For more detailed codes on MAGMA tool, please refer to [here](https://cloufield.github.io/GWASTutorial/09_Gene_based_analysis/)

### Example format

```
1) Single-cell data

The input format of single-cell data: Seurat-generated S4 object.

ctDRTF is fully compatiable with Seurat, a widely-used single-cell analysis tool.

2) MAGMA-result data (i.e., MAGMA_GWAS_data)

                SYMBOL      logP
1             HIST1H4L 37.543330
2                 DPYD 22.758030
3             HIST1H3I 21.643210
4               CACNB2 17.592524
5              CACNA1C 17.534141
6             PPP1R16B 17.196147
7             HIST1H4A 17.139482
8                 TCF4 16.818929
9               SFMBT1 16.348761
10               SFTA2 15.954599

````

