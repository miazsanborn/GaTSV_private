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
##y <- factor(df$sv_class, levels = c(0, 1), labels = c("germline", "somatic"))
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
  gamma = 0.1,
  probability = TRUE
)

# predicted probabilities on train set
probs_train <- attr(predict(svm_new, X_train, probability = TRUE), "probabilities")[, "somatic"]

# get optimal cutoff
cutoffs <- seq(0, 1, by = 0.001)
best_cutoff <- 0
best_score <- -Inf
ppv_values <- numeric(length(cutoffs))  

for (i in seq_along(cutoffs)) {
  c <- cutoffs[i]
  preds <- factor(ifelse(probs_train > c, "somatic", "germline"), levels = c("somatic", "germline"))
  
  ppv <- ppv_vec(y_train, preds, event_level = "first")
  tpr <- sens_vec(y_train, preds, event_level = "first")  
  #tpr <- sensitivity(preds, y_train)
  
  ppv_values[i] <- ppv
  score <- ppv + tpr

  if (!is.na(score) && score > best_score) {
    best_score <- score
    best_cutoff <- c
  }
}

cat("Optimal cutoff (train set):", best_cutoff, "\n")

# use cut-off on test set (get ppv w yardstick)
probs_test <- attr(predict(svm_new, X_test, probability = TRUE), "probabilities")[, "somatic"]
preds_test <- factor(ifelse(probs_test > best_cutoff, "somatic", "germline"), levels = c("somatic", "germline"))

# def truth and estimate
results_tbl <- tibble(truth = y_test, estimate = preds_test, prob = probs_test)


# generate numbers
tpr_test  <- sens_vec(results_tbl$truth, results_tbl$estimate, event_level = "first")
spec_test <- spec_vec(results_tbl$truth, results_tbl$estimate, event_level = "first")
ppv_test  <- ppv_vec(results_tbl$truth, results_tbl$estimate, event_level = "first")
auc_test  <- roc_auc_vec(results_tbl$truth, results_tbl$prob, event_level = "first")

# results
cat("Test AUC:", auc_test, "\n")
cat("Test PPV:", ppv_test, "\n")
cat("Test TPR:", tpr_test, "\n")
cat("Test Specificity:", spec_test, "\n")
cat("Optimal cutoff:", best_cutoff, "\n")
