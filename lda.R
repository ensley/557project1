library(dplyr)
library(MASS)

train <- sample(1:nrow(iris), nrow(iris)/2)
df <- iris[train, ]
classes <- levels(df$Species)
N <- nrow(df)
K <- length(classes)

# prior
pi <- as.vector(table(df$Species)/nrow(df))

# build sample means
means <- as.data.frame(matrix(NA, nrow = K, ncol = ncol(df)))
colnames(means) <- colnames(df)
for (i in 1:length(classes)) {
  means$Species[i] <- classes[i]
  df.sub <- df[df$Species == classes[i], ]
  for (j in 1:(ncol(df)-1)) {
    means[i,j] <- mean(df.sub[ ,j])
  }
}

# build sample covariance
sigma <- matrix(0, nrow=ncol(df)-1, ncol=ncol(df)-1)
for (i in 1:length(classes)) {
  df.sub <- df[df$Species == classes[i], ]
  for (j in 1:nrow(df.sub)) {
    vec <- as.numeric(df.sub[j,-5] - means[i,-5])
    sigma <- sigma + 1/(N-K) * (vec %*% t(vec))
  }
}

disc.func <- function(x, k, means, sigma, prior) {
  if (!(k %in% 1:K)) {
    stop('k out of bounds')
  }
  mu.k <- as.numeric(means[k,-5])
  as.numeric(t(x) %*% solve(sigma) %*% mu.k - 1/2 * t(mu.k) %*% solve(sigma) %*% mu.k + log(pi[k]))
}

G <- function(x, means, sigma, prior) {
  vals <- rep(-Inf, K)
  for (i in 1:K) {
    vals[i] <- disc.func(x, i, means, sigma, prior)
  }
  which.max(vals)
}

prediction <- rep(NA, nrow(iris))
for (i in 1:nrow(iris)) {
  prediction[i] <- G(as.numeric(iris[i,-5]), means, sigma, pi)
}
table(as.integer(iris$Species)[-train] == prediction[-train])

z <- lda(Species ~ ., iris, prior=pi, subset=train)
table(iris$Species[-train] == predict(z, iris[-train, ])$class)
