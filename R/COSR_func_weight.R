
## 2024-04-15
## Author: Yunlong Ma
## E-mail: glb-biotech@zju.edu.cn


#@' 3) Identifying cell type-specific regulons relevant to disease

#@' magma result processing
#@' magma_results <- magma_results %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
#@' MAGMA_GWAS_data <- magma_results[,c(10,11,8)]
#@' MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
#@' MAGMA_GWAS_data header: SYMBOL, logP, ZSTAT
#@' data_regulons1: TF-gene pairs matrix
#@' tf_left: all tf names
#@' MC_num: Set 1000 times of running MC_JSI_score_func()
#@' data_s1: matrix of genes and TFs specificity scores across cell types
#@' theta range from 0.1 ~ 1, default set to 0.5
#@' mi(power): expand the specificity difference between cell types, default set to 1
#@' mode: mo=1 indicates use magma-based z-scores as weights; mo=0 indicates no weights

COSR_func_weight <- function(tf_left=tf_left,
                      data_s1=data_s1,
                      data_regulons1=data_regulons1,
                      MAGMA_GWAS_data = MAGMA_GWAS_data,
                      MC_num = 1000,
                      Gene_num=500,
                      theta=0.5,
                      mi=1,
                      mo=1){
  
  #Name all regulons
  regulon_names <- data.frame(paste(tf_left,"_regulon",sep = ""))
  #tf_left: the vector of all TFs 
  #Specificity score of each regulon in a particular cell type
  #Output file
  #Collecting specificity scores
  Final_regulon_s <- data.frame(matrix(ncol = 1, nrow = length(tf_left)))
  names(Final_regulon_s)<-"ID_regulon"
  Final_regulon_s$ID_regulon <- tf_left
  
  #Collecting specificity*JSI for each regulon-disease link
  Final_regulon_ct_score<- data.frame(matrix(ncol = 1, nrow = length(tf_left)))
  names(Final_regulon_ct_score)<-"ID_regulon"
  Final_regulon_ct_score$ID_regulon <- tf_left
  
  #Collecting MC P value for each regulon-disease link
  Final_regulon_ct_mc_p<- data.frame(matrix(ncol = 1, nrow = length(tf_left)))
  names(Final_regulon_ct_mc_p)<-"ID_regulon"
  Final_regulon_ct_mc_p$ID_regulon <- tf_left
  
  #Obtain all names of the cell types
  all_celltype_names <- unique(data_s1[,3])
  
  
  #max-min scaling the specificity score

  for (m in all_celltype_names ){
    
    data_s1$scores[which(data_s1$celltypes == m)] <- max_min_scale(data_s1$scores[which(data_s1$celltypes == m)])
    
  }

  #Adjusted the other cell types
  #gene_list <- unique(data_s1$genes)
  
  #adj_score_all <- c()
  #for (x in 1:length(data_s1$genes)) {
  #temp_score <- data_s1[x,]$scores
  #sum_score <- sum(data_s1$scores[which(data_s1$genes==data_s1[x,]$genes)])
  #adj_score <- (temp_score+1)/(sum_score+1)
    
  #adj_score_all <- c(adj_score_all,adj_score)
  #}
  
  #collect all the adjusted specificity score
  #data_s1$adj_score <- adj_score_all
  

   "
   Using the COSR method to calculate the specificity scores (S) of TF-regulons;

   Using disease-specific genes and regulons weighted by S to identify cell type-specific TFs relevant to disease

   "
  #Open progress bar
  pb <- txtProgressBar(style=3)
  start_time <- Sys.time() ##record start time
     
  
  total_run <- length(all_celltype_names)*length(tf_left)
  count <- 0
  #COSR and JSI interaction analysis
  for (i in 1:length(all_celltype_names)){
    regulon_ct_score <- c()
    regulon_ct_mc_p <- c()
    regulon_s_all <-c()
 
    for (j in 1:length(tf_left)){
      
      #extracting the gene list of each regulon
      M1 <- data_regulons1[which(data_regulons1$regulons.meta.tf == tf_left[j]),]
      M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
      
      
      #all specificty score and z score of all regulon genes
      data_s1_sub <- data_s1[which(data_s1[,3] == all_celltype_names[i]),]
      
      #annotating magma z-score
      temp <- MAGMA_GWAS_data[which(MAGMA_GWAS_data$SYMBOL %in% data_s1_sub$genes),]
      temp <- temp[!duplicated(temp[,c("SYMBOL")]),]
      
      #overlap data_s1 regulon genes with magma genes
      data_s1_sub<- data_s1_sub[which(data_s1_sub$genes %in%temp$SYMBOL),]
      
      #match() function 
      #data_set for all regulon genes specificity and z scores
      data_s1_sub$magma_zscore <- temp$ZSTAT[match(data_s1_sub$genes, temp$SYMBOL )]
      
      
      #extracting the specificity score of each regulon
      data_s_M1 <- data_s1_sub[!is.na(match(data_s1_sub[,1], M1_regulon)),]
      #data_s_M1$anno <- rep("Gene",n_num,length(data_s_M1[,1]))
      data_s_M1$anno <- rep("Gene", length(data_s_M1[,1]))
      data_s_M1$anno[which(data_s_M1[,1] == M1_regulon[1])] <- "TF"
      
      #annotation MAGMA z-score
      #data_s_M1<- data_s_M1[which(data_s_M1$genes %in% MAGMA_GWAS_data$SYMBOL),]
      #data_s_M1$magma_zscore <- MAGMA_GWAS_data$ZSTAT[which(MAGMA_GWAS_data$SYMBOL %in% data_s_M1$genes)]
      
      ##choose analysis mode
      ni = mo
      
      #Calculating the module specificity score for TF in each regulon
      #tf_s_z <- data_s_M1[,c("adj_score","magma_zscore")][which(data_s_M1[,1] == M1_regulon[1]),]
      tf_s_z <- data_s_M1[,c("scores","magma_zscore")][which(data_s_M1[,1] == M1_regulon[1]),]
      tf_w <- as.numeric((tf_s_z[1])^mi*(tf_s_z[2])^ni)
      
      if(is.na(tf_w)){
        
        gene_s_z <- data_s_M1[,c("scores","magma_zscore")][which(data_s_M1[,1] != M1_regulon[1]),]
        gene_w <- as.numeric((gene_s_z[,1])^mi*(gene_s_z[,2])^ni)
        ave_s <- sum(gene_w)/(length(gene_w)+1)
        
        tf_w <- 0
        
        #theta = 0.5  #theta range from 0.1 ~ 1, default set to 0.5
        regulon_s <- as.numeric(tf_w) + as.numeric(theta*ave_s) #regulon-specific score for each cell type
        
        
      } else{
        
        #Calculating the module specificity score for genes in each regulon
        #gene_s_z <- data_s_M1[,c("adj_score","magma_zscore")][which(data_s_M1[,1] != M1_regulon[1]),]
        gene_s_z <- data_s_M1[,c("scores","magma_zscore")][which(data_s_M1[,1] != M1_regulon[1]),]
        gene_w <- as.numeric((gene_s_z[,1])^mi*(gene_s_z[,2])^ni)
        ave_s <- mean(gene_w)
        
        #theta = 0.5  #theta range from 0.1 ~ 1, default set to 0.5
        regulon_s <- as.numeric(tf_w) + as.numeric(theta*ave_s) #regulon-specific score for each cell type
        
      }
      
      #Sum
      regulon_s_all <- c(regulon_s_all,regulon_s)
      
      #Calculating the Jaccard Similarity Index (JSI)
      top_genes <- MAGMA_GWAS_data$SYMBOL[1:Gene_num]
      inter_genes1 <- length(intersect(top_genes,M1_regulon))
      union_genes1 <- length(union(top_genes,M1_regulon))
      JSI1 <- inter_genes1/union_genes1 # Jaccard similarity index
      
      #Interaction: specificity*JSI for each regulon-disease link
      ct_score <- regulon_s*JSI1
      
      #print(paste0("Regulon ",tf_left[j]," ctDRTF score is: ",ct_score, sep=""))
      
      
      #MC simulation for random specificity*JSI scores
      #Function: MC_JSI_score_func()
      #MC_num = 1000
      tf_left_1 <- tf_left[which(tf_left!=tf_left[j])] #removing the targeted TF as controls
      len_of_regulon = length(M1_regulon)
      all_genes <- MAGMA_GWAS_data$SYMBOL
      MC_results <- replicate(MC_num,MC_JSI_score_func_weight(data_s1_sub,
                                                       tf_left_1, 
                                                       len_of_regulon,
                                                       all_genes,
                                                       Gene_num,
                                                       mi,
                                                       mo))
      
      #Calculating the MC p-values
      MC_p <- (1+length(MC_results[MC_results>ct_score]))/(1+length(MC_results))
      

      #Running
      print(paste("Running the regulon of ",tf_left[j], " for the cell type of ",all_celltype_names[i],sep = ""))

      #Running percent:
      count=count+1
      completed_percent <- count/total_run
      print(sprintf('Completed percent: %1.2f%%',100*completed_percent))
      #Real-time progress bar
      #print(paste("Runing percent: ",percent((i+j)/(length(all_celltype_names)*length(tf_left))),sep = ""))
      setTxtProgressBar(pb,(count)/total_run)
      

      #Saving results
      regulon_ct_score <- c(regulon_ct_score,ct_score)
      regulon_ct_mc_p <- c(regulon_ct_mc_p,MC_p)
      

    }

    #Collecting specificity scores
    regulon_s_all<- as.data.frame(regulon_s_all)
    names(regulon_s_all) <- all_celltype_names[i]
    Final_regulon_s <- cbind(Final_regulon_s,regulon_s_all)
    
    #Collecting MC P values
    regulon_ct_mc_p<- as.data.frame(regulon_ct_mc_p)
    names(regulon_ct_mc_p) <- all_celltype_names[i]
    Final_regulon_ct_mc_p <- cbind(Final_regulon_ct_mc_p,regulon_ct_mc_p)
    
    #Normalization
  
    regulon_ct_score_norm <- (regulon_ct_score-mean(regulon_ct_score))/sd(regulon_ct_score)
    
    #Alternative normalized method
    #max-min normalization
    #regulon_ct_score_norm <- max_min_scale(regulon_ct_score)
    
    #Collecting specificity*JSI for each regulon-disease link
    regulon_ct_score_norm <- as.data.frame(regulon_ct_score_norm)
    names(regulon_ct_score_norm) <- all_celltype_names[i]
    Final_regulon_ct_score <- cbind(Final_regulon_ct_score,regulon_ct_score_norm)
    

  }
  
  ##Record end time
  end_time <- Sys.time() 
  
  #Close progress bar
  close(pb)
  
  #Calculating the running time
  run_time <- end_time - start_time
  print(paste("Running time: ",run_time),sep = "")
  
  #Outputs
  #out_results <- list(ctDRTF_score = Final_regulon_ct_score, 
  #                     MC_p = Final_regulon_ct_mc_p,
  #                      regulon_specificity_s=Final_regulon_s)
  #
  out_results <- list(ctDRTF_score = Final_regulon_ct_score, 
                      MC_p = Final_regulon_ct_mc_p)
  
  # 
 return(out_results)

}




