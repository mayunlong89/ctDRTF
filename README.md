# ctDRTF
A computational method for identifying cell type-specific transcription factor relevant to complex diseases


We assume that if GWAS-identified disease-specific genes are concordantly activated in a cell type-specific TF-regulon, then the TF is more likely to have a pivotal role in disease via the given cell type. Thus, we design a computational framework ctDRTF (cell type-specific Disease-Relevant Transcription Factor) for performing regulatory network-based inference of the associations between TF-related regulons and disease-specific gene sets in a context-specific manner. The input of scDRTF includes multimodal matrix of both scRNA-seq and scATAC-seq data, and GWAS summary data for a quantitative trait or disease (case-control study). 

![Workflow](https://github.com/mayunlong89/ctDRTF_analysis_codes/blob/main/figures_1/Figure%204.png)


# Step 1
## @ 1) Constructing global TF-gene regulatory network
## single_cell: the input single-cell data--Seurat object.
## n_genes: the minimum number of genes in a given regulon`

`GRN_func(single_cell = single_cell,n_genes=n_genes)`

# Step 2
## @ 2) Calculating specificity scores
## single_cell is the input single-cell data.
## num_genes <- length(rownames(single_cell))

`COSR_pre_func(single_cell=single_cell)`

# Step 3
## @ 3) Identifying cell type-specific regulons relevant to disease
## MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
## data_regulons1: TF-gene pairs matrix
## tf_left: all tf names
## MC_num: Set 100 times of running MC_JSI_score_func()
## data_s1: matrix of genes and TFs specificity scores across cell types
## theta range from 0 ~ 1, default set to 0.5

`COSR_func(tf_left=tf_left,
          data_s1=data_s1,
          data_regulons1=data_regulons1,
          MAGMA_GWAS_data = MAGMA_GWAS_data,
          MC_num = 100,
          Gene_num=500,
          theta=0.5)`













