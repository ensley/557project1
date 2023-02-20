mylogr <- function(y, x, maxit = 100, tol = 1e-5) {
  y <- as.numeric(y) - 1
  X <- as.matrix(cbind(1, x))      # add the intercept term to the predictors
  
  # constants
  N <- length(y)        # number of observations
  p <- ncol(X)          # number of predictors
  continue <- T
  i <- 1
  
  betas <- matrix(0, nrow = maxit, ncol = p)
  
  while (continue && i < maxit) {
    beta <- betas[i,]
    p <- exp(X %*% beta) / (1 + exp(X %*% beta))
    W <- diag(as.vector(p * (1 - p)))
    z <- X %*% beta + solve(W) %*% (y - p)
    betas[i + 1,] <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    if (all(abs(betas[i + 1,] - betas[i,]) / abs(betas[i,]) < tol)) {
      continue <- F
    }
    i <- i + 1
  }
  
  return(betas[i,])
  
}

mylogr.cv <- function(y, x, k, maxit = 100, tol = 1e-5) {
  cvs <- rep(0, k)
  
  # constants
  N <- length(y)        # number of observations
  
  # assign indices for which group each observation belongs to
  kappa <- sample(rep(1:k, length = N))
  
  for (i in 1:k) {
    y.tr <- y[kappa != i]     # training responses; all but ith group
    x.tr <- x[kappa != i,]    # training predictors; all but ith group
    y.va <- y[kappa == i]     # validation responses; ith group
    x.va <- x[kappa == i,]    # validation predictors; ith group
    
    beta <- mylogr(y.tr, x.tr, maxit, tol)
    
    X.va <- as.matrix(cbind(1, x.va))
    
    prediction <- rep(NA, length(y.va))
    for (j in 1:length(prediction)) {
      prediction[j] <- round(1 / (1 + exp(-t(beta) %*% X.va[j,])))
    }
    cvs[i] <- sum((prediction - (as.numeric(y.va) - 1)) ^ 2)
  }
  
  return(sum(cvs) / N)
  
}

# read data in and format it
data <- readRDS("data/input.rds")
system.time(cv.err.logr <- mylogr.cv(data$class, data[, -10], 10))
