# Abstract: Elastic logit (Ridge/Lasso regression) via coordinate decent
# Author: Quan-En (Taner)
# Final update: 2021-04-25

library(magrittr)
source("./utils.R")

elastic_logit <- function(X,y,labda=0.01,alpha=1,epsilon=1e-9,init_par=NULL,max_iter=50,return_history=FALSE){
  # Normalize input variables
  X <- apply(X, MARGIN=2, FUN=variable_norm)
  X <- cbind(1, X)
  n_var <- ncol(X)
  n_obs <- nrow(X)
  
  if (is.null(init_par)){
    init_par <- rep(0, n_var)
  }else if(length(init_par) != ncol(X)){
    print("init_par's length is not correct");break
  }
  
  gama <- labda * alpha
  
  update_beta <- init_par
  beta_hat <- init_par
  update_beta <- update_beta + 1
  
  update_history <- vector(mode = 'list')
  loop_num <- 0
  
  while(norm(matrix(update_beta - beta_hat)) > epsilon){
    update_beta <- beta_hat
    loop_num <- loop_num + 1
    if (loop_num > max_iter){cat('Over',loop_num,' iteration','\n') ; break}
    
    px <- sigmoid(X, beta_hat)
    px <- ifelse(1 - px < 1e-5, 1, ifelse(px < 1e-5, 0, px))
    px_weight <- ifelse(px %in% c(0,1), 1e-10, px * (1 - px))
    z <- ( X %*% beta_hat + ((y - px) / px_weight))
    W <- (t(px_weight) %*% X**2)[1,] / n_obs
    
    
    for (par_index in 1:length(init_par)){
      w <- W[par_index]
      y_hat <- X[,-par_index] %*% beta_hat[-par_index]
      update_value <- (t(px_weight) %*% ((z - y_hat) * X[,par_index]))[1] / n_obs
      
      if (par_index == 1){
        beta_hat[par_index] <- update_value / w
      }else{
        threshold_value <- ifelse(gama >= abs(update_value), 0, update_value + ifelse(update_value > 0, - gama, gama))
        beta_hat[par_index] <- threshold_value / (labda - gama +  w)
      }
    }
    beta_hat <- round(beta_hat, 6)
    if (return_history){
      update_history = rlist::list.append(update_history, beta_hat)
    }
    
    if (loop_num %% 5 == 0){print(paste0('Iteration : ', loop_num))}
  }
  print(paste0('Iteration : ', loop_num))
  beta_hat %<>% `names<-`(paste0('beta', seq(0, length(beta_hat) - 1)))
  if (return_history){
    return(update_history)
  }else{
    return(beta_hat)
  }
}