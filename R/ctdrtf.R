
#@ 2024-04-15 v1.2.0

####---------------------
#@ 2024-03-15 v1.1.1
#@ 2024-01-01 v1.1.0
#@ 2023-11-03 v1.0.0
####---------------------

#@ Author: Yunlong Ma
#@ E-mail: glb-biotech@zju.edu.cn
#@ E-mail: yunlong.ma@pennmedicine.upenn.edu

#main function
#single_cell: the input single-cell data.
#n_genes: the minimum number of genes in a given regulon 
#MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
#Gene_num: The number of disease-specific genes, default set to 500
#MC_num: Set 100 times of Monte Carlo simulation
#theta range from 0 ~ 1, default set to 0.5
#mi: expand the specificity difference between cell types,default set to 2
#mode=1 indicates use magma-based z-scores as weights; mode=0 indicates no weights

ctdrtf <- function(single_cell = single_cell,
                           MAGMA_GWAS_data = MAGMA_GWAS_data,
                           n_genes= 10,
                           Gene_num = 500,
                           MC_num = 1000,
                           theta=0.5,
                           mi=2,
                           mode=1){
   
  #1) Constructing global TF-gene regulatory network

  grn_outputs <- GRN_func(single_cell,n_genes = n_genes)
  
  #grn_outputs <- list(grn=data_regulons1,tf_names = tf_left)
  

  #2) Calculating specicity scores
  data_s1 <- COSR_pre_func(single_cell) 
  
  #3) Identifying cell type-specific regulons relevant to disease
  
    
  final_results <- COSR_func_weight(tf_left=grn_outputs$tf_names,
                                      data_s1=data_s1,
                                      data_regulons1=grn_outputs$grn,
                                      MAGMA_GWAS_data = MAGMA_GWAS_data)

  #4) Outputs
  return(final_results)
}
