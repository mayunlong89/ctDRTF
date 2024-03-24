# ctDRTF
A computational method for identifying cell type-specific transcription factor relevant to complex diseases


We assume that if GWAS-identified disease-specific genes are concordantly activated in a cell type-specific TF-regulon, then the TF is more likely to have a pivotal role in disease via the given cell type. Thus, we design a computational framework ctDRTF (cell type-specific Disease-Relevant Transcription Factor) for performing regulatory network-based inference of the associations between TF-related regulons and disease-specific gene sets in a context-specific manner. The input of scDRTF includes multimodal matrix of both scRNA-seq and scATAC-seq data, and GWAS summary data for a quantitative trait or disease (case-control study). 

![Workflow](https://github.com/mayunlong89/ctDRTF_analysis_codes/blob/main/figures_1/Figure%204.png)


#Step one
##@ 1) Constructing global TF-gene regulatory network
##single_cell: the input single-cell data--Seurat object.
##n_genes: the minimum number of genes in a given regulon`

`GRN_func(single_cell = single_cell,n_genes=n_genes)`
