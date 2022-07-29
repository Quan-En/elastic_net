# Abstract: Elastic net (Ridge/Lasso regression) via coordinate decent
# Author: Quan-En, Li
# Final update: 2021-04-25

library(magrittr)
source("./utils.R")

elastic_net <- function(X,y,labda=0.01,alpha=1,epsilon=1e-4,init_par=NULL,max_iter=100,return_history=F){
  
  # Normalize input variables
  X <- apply(X, MARGIN=2, FUN=variable_norm)
  X <- cbind(1, X)
  
  # Initialing parameter
  if (is.null(init_par)){
    init_par <- rep(0, ncol(X))
  }else if(length(init_par) != ncol(X)){
    print("init_par's length is not correct")
    break
  }
  
  sample_size <- length(y)
  
  gama <- labda * alpha
  
  beta_hat <- init_par
  
  update_beta <- 1
  
  update_history <- vector(mode = 'list')
  
  loop_num <- 0
  
  while(norm(matrix(update_beta - beta_hat)) > epsilon){
    
    update_beta <- beta_hat
    loop_num <- loop_num + 1
    
    if (loop_num > max_iter){cat('Over',loop_num,' interation','\n') ; break}
    
    for (par_index in 1:length(init_par)){
      
      if (par_index == 1){
        
        beta_hat[par_index] <- mean(y - X[,-par_index] %*% beta_hat[-par_index])
        
      }else{
        
        part1 <- X[, par_index] %*% y / sample_size
        part2 <- (apply(X, MARGIN = 2, FUN = `%*%`, X[,par_index]) %*% beta_hat) / sample_size
        
        z <- part1[1] - part2[1] + beta_hat[par_index]
        
        threshold_value <- ifelse(gama >= abs(z), 0, z + ifelse(z > 0, - gama, gama))
        
        beta_hat[par_index] <- threshold_value / (1 + labda - gama)
        
      }
    }
    
    if (loop_num %% 5 == 0){print(paste0('Interation : ',loop_num))}
    
    if (return_history){update_history <- rlist::list.append(update_history,beta_hat)}
    
  }
  
  print(paste0('Interation : ',loop_num))
  
  if (return_history){
    return(update_history)
  }else{
    beta_hat %<>% `names<-`(paste0('beta',1:length(beta_hat)))
    return(beta_hat)
  }
  
}