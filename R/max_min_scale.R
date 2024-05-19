
#@ 2024-04-08
#@ Author: Yunlong Ma
#@ E-mail: glb-biotech@zju.edu.cn

##-----------------------max-min scale--
max_min_scale <- function(x){
  
  for (i in 1: length(x)) {
    
    if(is.na(x[i])){
      
      x[i] <- mean(na.omit(x))
      
    } }
  
  scale_vec <- (x-min(x))/(max(x)-min(x))
  
  return(scale_vec)
  
}