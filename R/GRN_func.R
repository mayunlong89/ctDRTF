
#@ 2023-11-03
#@ Author: Yunlong Ma
#@ E-mail: glb-biotech@zju.edu.cn

#@ 1) Constructing global TF-gene regulatory network
#single_cell: the input single-cell data--Seurat object.
#n_genes: the minimum number of genes in a given regulon 

GRN_func<- function(single_cell = single_cell,n_genes=n_genes){
  # Select variable features
  single_cell <- Seurat::FindVariableFeatures(single_cell, assay='RNA')
  # Initiate GRN object and select candidate regions
  single_cell <- Pando::initiate_grn(single_cell)
  # Scan candidate regions for TF binding motifs
  single_cell <- Pando::find_motifs(
    single_cell,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  "
  Integrating single-cell transcriptomics and chromatin accessibility to construct GRN
   
  "
   
  # Infer gene regulatory network
  grn_object <- Pando::infer_grn(single_cell)
  
  # Find modules
  grn_object <- Pando::find_modules(grn_object)

  # Module extraction
  regulons <- Pando::NetworkModules(grn_object)
  
  # Extracting GRN
  data_regulons <- data.frame(regulons@meta$tf,regulons@meta$target)
  
  # Regulons with genes less than 10 were removed
  temp <- data.frame(table(data_regulons[,1]))
  temp <- temp[which(temp$Freq > n_genes-1),]
  
  # Number of TFs that construct number-matched regulons
  tf_left <- as.vector(temp$Var1)
  
  # All regulons remaining in the GRN
  data_regulons1 <- data_regulons[!is.na(match(data_regulons$regulons.meta.tf,tf_left)),]
  
  # Outputs
  grn_outputs <- list(grn=data_regulons1,tf_names = tf_left)
  
  return(grn_outputs)
}