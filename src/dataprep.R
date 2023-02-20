# read data in and format it
data <-
  read.csv('data/breast-cancer-wisconsin.data',
           header = F,
           na.strings = '?')
names(data) <- c(
  'id',
  'thickness',
  'size.unif',
  'shape.unif',
  'adhesion',
  'size.cell',
  'nuclei',
  'chromatin',
  'nucleoli',
  'mitoses',
  'class'
)
data$class <- as.factor(data$class)
levels(data$class) <- c('benign', 'malignant')
data <- data[-which(is.na(data$nuclei)), -1]  # leave ID column out and remove NAs

saveRDS(data, "data/input.rds")
