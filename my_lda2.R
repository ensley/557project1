# read data in and format it
data <- read.csv('breast-cancer-wisconsin.data', header = F, na.strings = '?')
names(data) <- c('id','thickness','size.unif','shape.unif','adhesion','size.cell','nuclei','chromatin','nucleoli','mitoses','class')
data$class <- as.factor(data$class)
levels(data$class) <- c('benign','malignant')

mylda <- function(y, x) {

  # constants
  N <- length(y)              # number of observations
  Nk <- as.vector(table(y))   # number of observations belonging to each class
  classes <- levels(y)        # vector containing the classes
  K <- length(classes)        # number of different classes
  p <- ncol(x)                # number of predictors
  
  # prior distribution
  prior <- Nk/N
  
  # sample means
  means <- as.data.frame(matrix(NA, nrow = K, ncol = p+1))
  colnames(means) <- c('class', colnames(x))
  for (i in 1:K) {
    means[i,1] <- classes[i]
    data.subset <- as.matrix(x[y == classes[i], ])
    means[i,-1] <- colMeans(data.subset)
  }
  
  # sample covariance
  sigma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:K) {
    data.subset <- as.matrix(x[y == classes[i], ])
    sigma <- sigma + cov(data.subset) * (Nk[i] - 1)
  }
  
  sigma <- sigma/(N-K)
  
  # discriminant function
  disc.func <- function(x, k, means, sigma, prior) {
    if (!(k %in% 1:K)) {
      stop('k out of bounds')
    }
    mu.k <- as.numeric(means[k,-1])
    eig <- eigen(sigma)
    v <- t(eig$vectors) %*% mu.k
    d.inv <- solve(diag(eig$values))
    as.numeric(t(x) %*% eig$vectors %*% d.inv %*% v - 1/2 * t(v) %*% d.inv %*% v + log(prior[k]))
  }
  
  # returns argmax of discriminant function
  G <- function(x, means, sigma, prior) {
    vals <- rep(-Inf, K)
    for (i in 1:K) {
      vals[i] <- disc.func(x, i, means, sigma, prior)
    }
    which.max(vals)
  }
  
  return(list(means=means,
              sigma=sigma,
              prior=prior,
              disc.func=disc.func,
              G=G,
              x=x,
              y=y))
}

data <- data[-which(is.na(data$nuclei)),-1] # leave ID column out and remove NAs
m <- mylda(data$class, data[ ,-10]) # about 95% correct classification


mylda.cv <- function(y, x, k) {
  cvs <- rep(0, k)

  # constants
  N <- length(y)              # number of observations
  Nk <- as.vector(table(y))   # number of observations belonging to each class
  classes <- levels(y)        # vector containing the classes
  K <- length(classes)        # number of different classes
  p <- ncol(x)                # number of predictors
  
  kappa <- sample(rep(1:k, length = N))
  
  for(i in 1:k) {
    y.tr <- y[kappa != i]
    x.tr <- x[kappa != i, ]
    y.va <- y[kappa == i]
    x.va <- x[kappa == i, ]
    m <- mylda(y.tr, x.tr)
    prediction <- rep(NA, length(y.va))
    for(j in 1:length(prediction)) {
      prediction[j] <- m$G(as.numeric(x.va[j, ]), m$means, m$sigma, m$prior)
    }
    cvs[i] <- sum((prediction - as.numeric(y.va))^2)
  }
  
  return(sum(cvs)/N)
}

system.time(cv.err.lda <- mylda.cv(data$class, data[ ,-10], 10))
