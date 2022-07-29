library(magrittr)


# Normalize
variable_norm <- function(x){return( (x - mean(x)) / sd(x) )}

# Sigmoid
sigmoid <- function(x,par_hat){
  lin_comb = (-1) * x %*% par_hat
  return(1 / (1 + exp(lin_comb)))
}

# Repeat experiment
repeat_experi <- function(x,y,FUN,experi_limit=50,repeat_cr=3,...){
  result_list <- vector(mode = 'list')
  stop_cr <- 0
  repeat{
    stop_cr = stop_cr + 1
    result_list <- rlist::list.append(.data = result_list,round(FUN(X,y,...),digits = 6))
    if (stop_cr > experi_limit){
      cat('Over',stop_cr,'time');break
    }else if (sum(duplicated(result_list)) > repeat_cr){
      cat('Number of experiments: ',stop_cr,'\n');break
    }
  }
  return(unique(result_list[duplicated(result_list)]))
}

# Calculate AIC and BIC
get_abic <- function(beta_coef,x,y,dist='gaussian'){
  X = apply(X = x,MARGIN = 2,variable_norm)
  X = cbind(1, X)
  n_obs = nrow(X)
  if (dist == 'gaussian'){
    pred = X %*% beta_coef
    rss = sum((pred - y)**2)
    k = sum(beta_coef != 0)
    model_aic = 2*k + n_obs*log(rss/n_obs)
    model_bic = (k+1)*log(n_obs) + n_obs*log(rss/n_obs)
  }else if (dist == 'binomial'){
    pred = sigmoid(X,beta_coef)
    logliklihood = sum(log(pred[y == 1]))+sum(log(1 - pred[y == 0]))
    k = sum(beta_coef != 0)
    model_aic = 2*k - 2*logliklihood
    model_bic = k*log(n_obs) - 2*logliklihood
  }
  return(list('AIC'=model_aic,'BIC'=model_bic))
}