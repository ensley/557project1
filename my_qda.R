# read data in and format it
data <- read.csv('breast-cancer-wisconsin.data', header = F, na.strings = '?')
names(data) <- c('id','thickness','size.unif','shape.unif','adhesion','size.cell','nuclei','chromatin','nucleoli','mitoses','class')
data$class <- as.factor(data$class)
levels(data$class) <- c('benign','malignant')

myqda <- function(data, resp.idx = 1) {
  train <- sample(1:nrow(data), floor(nrow(data)/2 + 1)) # split data set in half
  data.train <- data[train, ]
  data.test.df <- data[-train, ]
  data.test <- as.matrix(data.test.df[ ,-resp.idx])
  
  # constants
  N <- nrow(data.train)
  Nk <- as.vector(table(data.train[ ,resp.idx]))
  classes <- levels(data.train[ ,resp.idx])
  K <- length(classes)
  p <- ncol(data.train) - 1 # number of predictors
  
  # prior distribution
  prior <- Nk/N
  
  # sample means
  means <- as.data.frame(matrix(NA, nrow = K, ncol = ncol(data.train)))
  colnames(means) <- colnames(data.train)
  for (i in 1:K) {
    means[i,10] <- classes[i]
    data.subset <- as.matrix(data.train[data.train$class == classes[i],-resp.idx])
    means[i,-10] <- colMeans(data.subset)
  }
  
  # sample covariance
  sigma <- array(0, dim = c(p, p, 2))
  for (i in 1:K) {
    data.subset <- data.train[data.train[ ,resp.idx] == classes[i], ]
    data.mat <- as.matrix(data.subset[ ,-resp.idx])
    sigma[ , ,i] <- cov(data.mat)
  }
  
  # discriminant function
  disc.func <- function(x, k, means, sigma, prior) {
    if (!(k %in% 1:K)) {
      stop('k out of bounds')
    }
    sig <- sigma[ , ,k]
    mu.k <- as.numeric(means[k,-resp.idx])
    eig <- eigen(sig)
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
    classes[which.max(vals)]
  }
  
  # prediction on test data
  prediction <- rep(NA, nrow(data.test))
  for (i in 1:length(prediction)) {
    prediction[i] <- G(data.test[i, ], means, sigma, prior)
  }
  table(prediction == data.test.df$class)/sum(table(prediction == data.test.df[ ,resp.idx]))
}

data <- data[-which(is.na(data$nuclei)),-1] # leave ID column out and remove NAs
myqda(data, 10) # about 95% correct classification
