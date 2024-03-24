
#@ 2023-11-03
#@ Author: Yunlong Ma
#@ E-mail: glb-biotech@zju.edu.cn


#@ 2) Calculating specificity scores
#single_cell is the input single-cell data.
#num_genes <- length(rownames(single_cell))

COSR_pre_func <- function(single_cell = single_cell){

  #Calculating the specificity scores of TFs and targeting genes across all cell types
  
  #All genes in single-cell data 
  num_genes <- length(rownames(single_cell))

  #COSG::cosg()
  celltype_markers <- COSG::cosg(single_cell,
                                 groups = "all",
                                 assay = "RNA",
                                 slot = "data",
                                 mu=1,
                                 n_genes_user = num_genes) #number of genes used for analysis
  
  

  #extracting cell type names
  all_celltype_names <- colnames(celltype_markers$names)
  
  ##-----Extracting the specificity score of each gene or TF in a given regulon
  #Empty data.frame
  data_s1 <- data.frame(matrix(ncol = 3, nrow = 0))
  names(data_s1) <- c("genes","scores","celltypes")
  for (i in 1:length(all_celltype_names)){
    
    data_s <- data.frame(celltype_markers$names[i],celltype_markers$scores[i])
    data_s[,3] <- rep(all_celltype_names[i],length(data_s[,1]))
    names(data_s) <- c("genes","scores","celltypes")
    
    #collecting all the cell types
    data_s1 <- rbind(data_s1,data_s)
  }
  
  #If specificity score = -1, transforming -1 to 0
  data_s1[,2][which(data_s1[,2] == -1)] <- 0
  
  return(data_s1)

}
#@ 2023-11-03
#@ Author: Yunlong Ma
#@ E-mail: glb-biotech@zju.edu.cn


#@ 2) Calculating specificity scores
#single_cell is the input single-cell data.
#num_genes <- length(rownames(single_cell))

COSR_pre_func <- function(single_cell = single_cell){

  #Calculating the specificity scores of TFs and targeting genes across all cell types
  
  #All genes in single-cell data 
  num_genes <- length(rownames(single_cell))

  #COSG::cosg()
  celltype_markers <- COSG::cosg(single_cell,
                                 groups = "all",
                                 assay = "RNA",
                                 slot = "data",
                                 mu=1,
                                 n_genes_user = num_genes) #number of genes used for analysis
  
  

  #extracting cell type names
  all_celltype_names <- colnames(celltype_markers$names)
  
  ##-----Extracting the specificity score of each gene or TF in a given regulon
  #Empty data.frame
  data_s1 <- data.frame(matrix(ncol = 3, nrow = 0))
  names(data_s1) <- c("genes","scores","celltypes")
  for (i in 1:length(all_celltype_names)){
    
    data_s <- data.frame(celltype_markers$names[i],celltype_markers$scores[i])
    data_s[,3] <- rep(all_celltype_names[i],length(data_s[,1]))
    names(data_s) <- c("genes","scores","celltypes")
    
    #collecting all the cell types
    data_s1 <- rbind(data_s1,data_s)
  }
  
  #If specificity score = -1, transforming -1 to 0
  data_s1[,2][which(data_s1[,2] == -1)] <- 0
  
  return(data_s1)

}