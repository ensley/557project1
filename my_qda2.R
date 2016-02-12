# read data in and format it
data <- read.csv('breast-cancer-wisconsin.data', header = F, na.strings = '?')
names(data) <- c('id','thickness','size.unif','shape.unif','adhesion','size.cell','nuclei','chromatin','nucleoli','mitoses','class')
data$class <- as.factor(data$class)
levels(data$class) <- c('benign','malignant')
data <- data[-which(is.na(data$nuclei)),-1] # leave ID column out and remove NAs

myqda <- function(y, x) {
  
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
  sigma <- array(0, dim = c(p, p, 2))
  for (i in 1:K) {
    data.subset <- as.matrix(x[y == classes[i], ])
    sigma[ , ,i] <- cov(data.subset)
  }

  # discriminant function
  disc.func <- function(x, k, means, sigmas, prior) {
    if (!(k %in% 1:K)) {
      stop('k out of bounds')
    }
    sigma <- sigmas[ , ,k]
    mu.k <- as.numeric(means[k,-1])
    eig <- eigen(sigma)
    v <- t(eig$vectors) %*% (x - mu.k)
    d.inv <- solve(diag(eig$values))
    as.numeric(-1/2 * sum(log(eig$values)) - 1/2 * t(v) %*% d.inv %*% v + log(prior[k]))
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


m <- myqda(data$class, data[ ,-10]) # about 95% correct classification


myqda.cv <- function(y, x, k) {
  cvs <- rep(0, k)
  
  # constants
  N <- length(y)              # number of observations
  Nk <- as.vector(table(y))   # number of observations belonging to each class
  classes <- levels(y)        # vector containing the classes
  K <- length(classes)        # number of different classes
  p <- ncol(x)                # number of predictors
  
  # assign indices for which group each observation belongs to
  kappa <- sample(rep(1:k, length = N))
  
  # leave one different group out of each fit
  for(i in 1:k) {
    y.tr <- y[kappa != i]     # training responses; all but ith group
    x.tr <- x[kappa != i, ]   # training predictors; all but ith group
    y.va <- y[kappa == i]     # validation responses; ith group
    x.va <- x[kappa == i, ]   # validation predictors; ith group
    # fit the QDA model to the training data
    m <- myqda(y.tr, x.tr)
    prediction <- rep(NA, length(y.va))
    for(j in 1:length(prediction)) {
      # predict the jth response using the fitted model and the validation predictors
      prediction[j] <- m$G(as.numeric(x.va[j, ]), m$means, m$sigma, m$prior)
    }
    # store squared error for this fit
    cvs[i] <- sum((prediction - as.numeric(y.va))^2)
  }
  # calculate mean squared prediction error
  return(sum(cvs)/N)
}

system.time(cv.err.qda <- myqda.cv(data$class, data[ ,-10], 10))
# looks like LDA does better than QDA (~0.3 vs ~0.4)
