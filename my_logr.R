# read data in and format it
data <- read.csv('breast-cancer-wisconsin.data', header = F, na.strings = '?')
names(data) <- c('id','thickness','size.unif','shape.unif','adhesion','size.cell','nuclei','chromatin','nucleoli','mitoses','class')
data$class <- as.factor(data$class)
levels(data$class) <- c('benign','malignant')
data <- data[-which(is.na(data$nuclei)),-1] # leave ID column out and remove NAs

train <- sample(1:683, 342)
data.train <- data[train, ]
data.test <- data[-train, ]

y <- as.numeric(data.train[ ,10]) - 1
X <- cbind(rep(1, length(y)), as.matrix(data.train[ ,-10]))
beta0 <- rep(0, 10)

maxit <- 100
betas <- matrix(0, nrow = maxit, ncol = 10)
continue <- T
iters <- 1
tol <- 1e-5


while(continue && iters < maxit) {
  beta <- betas[iters, ]
  p <- exp(X %*% beta)/(1 + exp(X %*% beta))
  W <- diag(as.vector(p * (1-p)))
  z <- X %*% beta + solve(W) %*% (y-p)
  betas[iters + 1, ] <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
  if(all(abs(betas[iters + 1, ] - betas[iters, ])/abs(betas[iters, ]) < tol)) {
    continue = F
  }
  iters = iters + 1
}

beta <- betas[iters, ]

glm(class ~ ., data = data, family = binomial, subset = train)
