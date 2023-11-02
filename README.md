# ctDRTF

We assume that if GWAS-identified disease-specific genes are concordantly activated in a cell type-specific TF-regulon, then the TF is more likely to have a pivotal role in disease via the given cell type. Thus, we design a computational framework ctDRTF (cell type-specific Disease-Relevant Transcription Factor) for performing regulatory network-based inference of the associations between TF-related regulons and disease-specific gene sets in a context-specific manner. The input of scDRTF includes multimodal matrix of both scRNA-seq and scATAC-seq data, and GWAS summary data for a quantitative trait or disease (case-control study). 

![Figure 1. Framework of ctDRTF](https://github.com/mayunlong89/ctDRTF/blob/main/figures_1/Figure%204.png)

This computational framework has been done and applied in the bioRvix 2023.

# scHOB Database
scHOB (single-cell Human Organoid Bank) is a multi-omic single-cell database that consists of both scRNA-seq and scATAC-seq data on 10 kinds of widely-adopted human organoids (i.e., brain, lung, heart, eye, liver & bile duct, pancreas, intestine, kidney, and skin) spanning 1,366,380 cells with 67 main cell types in 385 samples across 83 distinct protocols. see github [Link](https://github.com/mayunlong89/scHOB/tree/main)  see Website [link](https://schob.su-lab.org/)

# Citation
Ma et al. Integrated analysis of single-cell transcriptomic and epigenomic atlas of organoids uncovers core regulomes of human diseases. bioRvix, 2023
