# read data in and format it
data <- read.csv('breast-cancer-wisconsin.data', header = F, na.strings = '?')
names(data) <- c('id','thickness','size.unif','shape.unif','adhesion','size.cell','nuclei','chromatin','nucleoli','mitoses','class')
data$class <- as.factor(data$class)
levels(data$class) <- c('benign','malignant')
data <- data[-which(is.na(data$nuclei)),-1] # leave ID column out and remove NAs

X <- as.matrix(data[ ,-10])

X <- apply(X, 2, function(x) (x - mean(x))/sd(x))

eig <- eigen(cov(X))
eig$values/sum(eig$values)

V <- eig$vectors[ ,1:2]
newX <- as.matrix(X %*% V)

# mylda(data$class, newX)
mylda.cv(data$class, newX, 10)
m <- mylda(data$class, newX)

plot(newX, col=data$class, pch=20)

a0 <- log(m$prior[1]/m$prior[2]) - 1/2*as.matrix(m$means[1,2:3] + m$means[2,2:3]) %*% solve(m$sigma) %*% t(as.matrix(m$means[1,2:3] - m$means[2,2:3]))
a1 <- solve(m$sigma) %*% t(as.matrix(m$means[1,2:3] - m$means[2,2:3]))

abline(-a0/a1[2], -a1[1]/a1[2])


myqda.cv(data$class, newX, 10)
mylogr.cv(data$class, newX, 10)
