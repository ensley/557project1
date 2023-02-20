# read data in and format it
data <- readRDS("data/input.rds")

X <- as.matrix(data[,-10])
# center and scale the data
X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))

eig <- eigen(cov(X))
eig$values / sum(eig$values)

V <- eig$vectors[, 1:2]
newX <- as.matrix(X %*% V)

# calculate decision boundary line coefficients
m <- mylda(data$class, newX)
a0 <- log(m$prior[1] / m$prior[2]) - 1 / 2 * as.matrix(m$means[1, 2:3] + m$means[2, 2:3]) %*% solve(m$sigma) %*% t(as.matrix(m$means[1, 2:3] - m$means[2, 2:3]))
a1 <- solve(m$sigma) %*% t(as.matrix(m$means[1, 2:3] - m$means[2, 2:3]))

pca.lda <- mylda.cv(data$class, newX, 10)
pca.qda <- myqda.cv(data$class, newX, 10)
pca.logr <- mylogr.cv(data$class, newX, 10)

c(LDA = pca.lda, QDA = pca.qda, log.r = pca.logr)
