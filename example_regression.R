

# Import packages and model
library(magrittr)
source("./model/elastic_net.R")

# Load data
train_data <- read.csv("./data/train.csv")
test_data <- read.csv("./data/test.csv")

train_x <- train_data[,1:121]
train_y <- ifelse(train_data[['classe']] == "b'A'", 1, 0)

test_x <- test_data[,1:121]
test_y <- ifelse(test_data[['classe']] == "b'A'", 1, 0)

# Lasso regression
## Set alpha=1 <=> Lasso, alpha=0 <=> Ridge
example_result <- elastic_net(X=train_x, y=(2 * train_y) - 1, labda=0.01, alpha=1, epsilon=0.0001)


# Fine tune lambda

## Create empty list to store estimation result
all_beta_coef <- vector(mode="list")
## Create empty list to store model AIC/BIC
all_abic <- vector(mode="list")

candidate_lambda_seq <- seq(0.010, 0.020, 0.001)

## Initialize parameter
estimation_result <- rep(0, ncol(train_x) + 1)

for (candidate_lambda in candidate_lambda_seq){
  
  ## Build elastic regression iteratively and estimate parameter
  estimation_result <- elastic_net(
    X=train_x,
    y=(2 * train_y) - 1,
    labda=candidate_lambda,
    init_par=estimation_result,
    epsilon=0.0001
  )
  
  ## Calculate AIC/BIC
  abic <- get_abic(beta_coef=estimation_result, x=train_x, y=(2 * train_y) - 1, dist="gaussian")
  
  ## Store result
  all_beta_coef <- rlist::list.append(all_beta_coef, estimation_result)
  all_abic <- rlist::list.append(all_abic, abic)
  
  
  print(candidate_lambda)
  print(abic)
}

# Prediction

## Estimation
final_estimation <- elastic_net(
  X=train_x,
  y=(2 * train_y) - 1,
  labda=0.02,
  alpha=1,
  init_par=estimation_result,
  epsilon=0.0001
)

## predict
lin_pred <- ifelse(cbind(1, apply(X=test_x, MARGIN=2, FUN=variable_norm)) %*% final_estimation > 0, 1, -1)

## calculate accuracy
classification_table <- table(lin_pred, (2 * test_y) - 1)
acc <- sum(diag(classification_table)) / sum(classification_table)
print(acc)

# Visualization of estimation result

## collect estimation history
beta_history <- elastic_net(
  X=train_x,
  y=(2 * train_y) - 1,
  labda=0.02,
  alpha=1,
  epsilon=0.0001,
  return_history=TRUE
)

beta_history <- do.call(rbind, beta_history)
beta_history <- rbind(rep(0, ncol(beta_history)), beta_history)

## plot
matplot(
  beta_history,
  pch=1,
  type='l',
  col=4,
  xlab='Iteration',
  ylab='Coefficient',
  lty=1,
  main='Linear regression coefficient (Lasso)'
)

# Compare to package's result
pkg_result <- glmnet::glmnet(x=model.matrix(train_y ~ ., data=data.frame(train_y, apply(train_x, 2, variable_norm))),
                             y=train_y, family='gaussian', lambda=0.02, intercept=TRUE, alpha=1)
pkg_result[['beta']]
