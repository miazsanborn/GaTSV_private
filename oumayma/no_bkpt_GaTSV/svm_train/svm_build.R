.libPaths(c("/data/beroukhim1/oumayma/Rlibs", .libPaths()))

# libraries
library(e1071)
library(caret)
library(pROC)
library(data.table)
library(tibble)
library(yardstick)
set.seed(123)

# tcga test-train dataset
df <- readRDS("/data/beroukhim1/oumayma/no_bkpt_GaTSV/svm_train/filtered_TCGA_scaled.rds")
df <- as.data.table(df)

# remove breakpoint features
exclude_features <- c("log_homlen", "log_insertion_len", "hom_gc", "insertion_gc", "sv_class")
keep_cols <- setdiff(names(df), exclude_features)

X <- as.matrix(df[, ..keep_cols])
colnames(X) <- keep_cols

# labels (0 = germline, 1 = somatic)
y <- factor(df$sv_class, levels = c(1, 0), labels = c("somatic", "germline"))


# split dataset: 2 train : 1 test
trainIndex <- createDataPartition(y, p = 0.66, list = FALSE)
X_train <- X[trainIndex, ]
y_train <- y[trainIndex]
X_test  <- X[-trainIndex, ]
y_test  <- y[-trainIndex]

# use gatsv params
svm_new <- svm(
  x = X_train,
  y = y_train,
  kernel = "radial",
  cost = 10,
  gamma = 0.05,
  probability = TRUE
)

save(svm_new, X_train, y_train, X_test, y_test, file = "/data/beroukhim1/oumayma/no_bkpt_GaTSV/GaTS_BpA_031826.rda")
