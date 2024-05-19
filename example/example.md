
# How to run example using ctDRTF

```r
#Single-cell multi-omic data on pbmc_10x (Illumina Official Website)
single_cell <- readRDS("pbmc_10x_example.rds")

#MAGMA-based gene-based analysis of GWAS on monocyte count (UKBiobank ID:ieu-b-31)
magma_example <- readRDS("magma_example.rds")

#running ctDRTF
data_example <- ctdrtf(single_cell = single_cell,
                               MAGMA_GWAS_data = magma_example,
                               n_genes = 10,
                               Gene_num = 500,
                               MC_num=1000,
                               theta=0.5)

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
magma_results <- read.table("magma.genes.out",header = TRUE)
magma_results <- magma_results %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
MAGMA_GWAS_data <- magma_results[,c(10,11,8)]

##---magma result processing
#MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
#header of MAGMA_GWAS_data: SYMBOL, logP, ZSTAT

```
#For more detailed codes on MAGMA tool, please refer to [here](https://cloufield.github.io/GWASTutorial/09_Gene_based_analysis/)

### Example format

```
1) Single-cell data

The input format of single-cell data: Seurat-generated S4 object.

ctDRTF is fully compatiable with Seurat, a widely-used single-cell analysis tool.

2) MAGMA-result data (i.e., MAGMA_GWAS_data)

              SYMBOL      logP  ZSTAT
1             PON1 10.899319 6.6721
2             PON2  8.570943 5.8352
3             PON3  4.386613 3.9381
4         CDKN2AIP  4.130850 3.7944
5             CASR  4.083372 3.7672
6             ASB4  3.919446 3.6719
7          UBASH3A  3.786084 3.5927
8             FUT6  3.776270 3.5868
9             WBP2  3.705379 3.5440
10         CCDC170  3.616292 3.4895
11           SFTPD  3.550830 3.4490
12            DKK1  3.442132 3.3809
13           LYRM7  3.429761 3.3730
14           FRG1B  3.384881 3.3445
15           NTNG1  3.373680 3.3373
16         TRAPPC9  3.360842 3.3291
17            NFIB  3.280420 3.2772
18           NUDT6  3.264058 3.2665
19         ADIPOR1  3.237014 3.2488
20           RPP30  3.229619 3.2440
21           DOCK9  3.225921 3.2416
22           ZFPM2  3.222283 3.2392
23           CDH24  3.207741 3.2296
24           FNIP2  3.187053 3.2160
25          CAPZA1  3.181358 3.2122
26            ST7L  3.135240 3.1816
27           PRR15  3.132144 3.1795
28          ANKRD1  3.112433 3.1663
29         AKIRIN1  3.065779 3.1349
30           DMXL2  3.054748 3.1275
31          PSMB11  3.016040 3.1012
32    RP11-527L4.2  3.000000 3.0902
33           ICAM2  2.985102 3.0800
34           STX19  2.968713 3.0688
35          SH3GL3  2.956481 3.0604
36         SLC35C1  2.954521 3.0590
37          CAMK2D  2.949002 3.0552
38          BTBD16  2.922814 3.0371
39          DZIP1L  2.917322 3.0332
40           PANK1  2.914068 3.0310
41           ZNF66  2.912432 3.0299
42        SIGLECL1  2.889141 3.0136
43        GTF2IRD1  2.886725 3.0119
44           PROS1  2.881934 3.0086
45          ZNF609  2.871375 3.0012
46           AP1M2  2.867100 2.9982
47           BSPRY  2.853221 2.9884
48            PAX1  2.850842 2.9868
49          CYB5R1  2.847253 2.9842
50           TMCO3  2.838632 2.9781
51           LARP6  2.836660 2.9768
52           NR6A1  2.834874 2.9755
53           SPRR4  2.833363 2.9744
54      AL162389.1  2.830002 2.9720
55          COX4I2  2.825068 2.9686
56           NSUN3  2.814373 2.9610
57          TRIM61  2.804266 2.9538
58          SPATA5  2.794471 2.9468
59           MARK4  2.782306 2.9382
60          LGALS2  2.775415 2.9332
61             MYC  2.756714 2.9198
62         C1orf64  2.742609 2.9097
63          SPRR2F  2.739166 2.9072
64            PGA3  2.729484 2.9003
65         C4orf22  2.704719 2.8823
66         CACNA1H  2.690604 2.8721
67        C17orf72  2.675327 2.8610
68           NR1H4  2.674361 2.8602
69            CINP  2.673050 2.8593
70           CLVS1  2.661962 2.8512
71            GPX1  2.655156 2.8462
72          UNC13D  2.645047 2.8388
73         DCUN1D2  2.644970 2.8387
74            SKA1  2.642370 2.8368
75           NUDT4  2.626334 2.8250
76          ARL13B  2.620821 2.8209
77          CDKN2D  2.617425 2.8184
78           MS4A5  2.607497 2.8111
79          SFTPA1  2.604464 2.8088
80           TRIP4  2.592813 2.8002
81   RP11-451M19.3  2.592099 2.7996
82           KLRK1  2.586231 2.7953
83          TRIM65  2.585444 2.7947
84         C1orf65  2.583593 2.7933
85           ZNF85  2.562139 2.7773
86             AMT  2.560383 2.7760
87        HSD17B12  2.551603 2.7694
88           MAPK6  2.550907 2.7689
89           TPRX1  2.521030 2.7464
90         SLC46A1  2.519145 2.7450
91          SNAP25  2.504442 2.7338
92            FGF2  2.496005 2.7274
93          CYP4B1  2.488023 2.7214
94           PNLIP  2.485466 2.7194
95            EEF2  2.480579 2.7157
96         ADCYAP1  2.473221 2.7101
97          ELOVL1  2.471148 2.7085
98            PUM2  2.467424 2.7056
99            GKN2  2.465860 2.7044
100       KIAA1009  2.460911 2.7007

````

