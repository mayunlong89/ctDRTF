
#@ 2023-11-03
#@ Author: Yunlong Ma
#@ E-mail: glb-biotech@zju.edu.cn


#@ 4) Calculating the P value for each regulon-disease link
#MC Function
#data_s1_sub: The input of matrix containing the specificity score of genes in each cell type
#tf_left: All the high-qualified TFs that passed the quality control
#tf_left_1 <- tf_left[which(tf_left!=tf_left[j])] #removing the targeted TF as controls
#len_of_regulon = length(M1_regulon), the number of genes and TFs in a given regulon
#all_genes: All genes from phenotype-associations based on MAGMA,ie.,length(MAGMA_GWAS_data$SYMBOL)
#Gene_num: The number of disease-specific genes, default set to 500

#@ Random-based specificity*JSI score for each regulon
MC_JSI_score_func_weight<- function(data_s1_sub = data_s1_sub,
                             tf_left_1 = tf_left_1,
                             len_of_regulon = len_of_regulon,
                             all_genes = all_genes, 
                             Gene_num = 500,
                             theta=0.5){
  
  #Selecting one TF
  select_tf <- sample(tf_left_1,1)
  
  #Random selecting genes
  data_s1_sub2 <- data_s1_sub[which(data_s1_sub$genes!=select_tf),]
  all_id <- rownames(data_s1_sub2)
  len_of_regulon_removing_one <- len_of_regulon-1 #removing the one count of TF
  order_id <- sample(all_id,len_of_regulon_removing_one)
  sample_data <- data_s1_sub2[order_id,]
  
  
  #Curating one selected TF and genes for a random regulon
  regulon_genes_sample <- c(select_tf,sample_data$genes)
  
  #Calculating the TF and targeting genes specificity score
  tf_s_sample <- data_s1_sub[,c("scores","magma_zscore")][which(data_s1_sub$genes==select_tf),]
  tf_s_sample_w <- as.numeric(tf_s_sample[1]*tf_s_sample[2])
   
  #Calculating the module specificity score for genes in each regulon
  ave_s_sample <- mean(as.numeric(sample_data[,2]*sample_data[,4]))
  
  #Calculating the module specificity score
  regulon_s_sample <- tf_s_sample_w + theta*ave_s_sample
  
  #Selecting the number-matched genes from background genes
  bg_genes <- sample(all_genes, Gene_num)
  
  #Calculating the Jaccard Similarity Index(JSI)
  inter_genes <- length(intersect(bg_genes,regulon_genes_sample))
  union_genes <- length(union(bg_genes,regulon_genes_sample))  
  JSI_sample <- inter_genes/union_genes  
  
  #Interaction: specificity*JSI for each regulon-disease link
  ct_score_sample <- regulon_s_sample*JSI_sample 
  
  return(ct_score_sample)  
  
}