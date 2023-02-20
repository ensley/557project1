library(ggplot2)
library(ggthemes)

source("src/lda.R")
source("src/qda.R")
source("src/logr.R")
source("src/pca.R")

pc <- data.frame(x = 1:9, y = eig$values / sum(eig$values))

ggplot(pc, aes(x, y)) + 
  geom_point(size = 3, shape = 15) + 
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = 1:9) +
  labs(title = 'Principal Components', x = '', y = 'PC') +
  theme_few() + scale_color_few()
ggsave("pca.png", path="images", width=963, height=461, units="px", dpi=100)

pca.lda.data <- data.frame(x = newX[ ,1], y = newX[ ,2], class = data$class)

ggplot(pca.lda.data, aes(x, y)) +
  geom_point(aes(color = class), size = 3, alpha = 0.9) +
  geom_abline(intercept = -a0 / a1[2], slope = -a1[1] / a1[2]) +
  labs(title = 'Decision Boundary for LDA') +
  theme_few() + scale_color_solarized()
ggsave("lda.png", path="images", width=963, height=461, units="px", dpi=100)

cv.data <- data.frame(type = rep(c('LDA','QDA','LogReg'), each = 5), err = rep(0, 15), folds = rep(1:5*2, 3))
cv.pca.data <- data.frame(type = rep(c('LDA','QDA','LogReg'), each = 5), err = rep(0, 15), folds = rep(1:5*2, 3))

set.seed(2015)
for (i in 1:15) {
  if (cv.data$type[i] == "LDA") {
    cv.data$err[i] <- mylda.cv(data$class, data[, -10], cv.data$folds[i])
  } else if (cv.data$type[i] == "QDA") {
    cv.data$err[i] <- myqda.cv(data$class, data[, -10], cv.data$folds[i])
  } else if (cv.data$type[i] == "LogReg") {
    cv.data$err[i] <- mylogr.cv(data$class, data[, -10], cv.data$folds[i])
  }
}

for (i in 1:15) {
  if (cv.pca.data$type[i] == "LDA") {
    cv.pca.data$err[i] <- mylda.cv(pca.lda.data$class, pca.lda.data[, -3], cv.pca.data$folds[i])
  } else if (cv.pca.data$type[i] == "QDA") {
    cv.pca.data$err[i] <- myqda.cv(pca.lda.data$class, pca.lda.data[, -3], cv.pca.data$folds[i])
  } else if (cv.pca.data$type[i] == "LogReg") {
    cv.pca.data$err[i] <- mylogr.cv(pca.lda.data$class, pca.lda.data[, -3], cv.pca.data$folds[i])
  }
}


ggplot(cv.data, aes(folds, err)) +
  geom_line(aes(color = type), size = 1) +
  geom_point(aes(color = type), size = 3) +
  ylim(0, 0.055) + labs(title = 'Cross Validation, Original Data') +
  theme_few() + scale_color_solarized()
ggsave("cv-nopc.png", path="images", width=963, height=461, units="px", dpi=100)


ggplot(cv.pca.data, aes(folds, err)) +
  geom_line(aes(color = type), size = 1) +
  geom_point(aes(color = type), size = 3) +
  ylim(0, 0.055) + labs(title = 'Cross Validation, 2 PCs') +
  theme_few() + scale_color_solarized()
ggsave("cv-pc.png", path="images", width=963, height=461, units="px", dpi=100)
