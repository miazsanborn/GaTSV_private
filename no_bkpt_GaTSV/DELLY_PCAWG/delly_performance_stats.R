# libraries
library(e1071)
library(caret)
library(pROC)
library(data.table)
library(tibble)
library(yardstick)
library(dplyr)

# replicates
n_reps <- 10
set.seed(1234)
seeds <- sample(1:10000, n_reps, replace = FALSE)

# 1. Load the unscaled combined data
unscaled_data <- readRDS("/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/processed/combined_delly.rds")
unscaled_data <- as.data.table(unscaled_data)

# 2. Apply the exact same filter used during scaling to drop the extra rows
unscaled_data <- unscaled_data[SPAN >= 1e3 | SPAN == -1, ]

# 3. Load your scaled data
df <- readRDS("/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/delly_test_scaled.rds")
df <- as.data.frame(df) # Safe to use data.frame here!

# THE FIX: Replace all NA values with 0 so the SVM doesn't drop rows and crash
df[is.na(df)] <- 0

# Sanity check: Ensure they are now the exact same length before we blindly paste
if (nrow(df) != nrow(unscaled_data)) {
  stop(paste("Row mismatch! Scaled data has", nrow(df), "rows, but filtered combined data has", nrow(unscaled_data)))
}

# 4. Create the sv_class column mapping (SOMATIC = 1, GERMLINE = 0)
df$sv_class <- ifelse(unscaled_data$truth == "SOMATIC", 1, 0)

# Now proceed exactly as before...
exclude_features <- c("log_homlen", "log_insertion_len", "hom_gc", "insertion_gc", "sv_class")
keep_cols <- setdiff(names(df), exclude_features)
X <- df[, keep_cols]
colnames(X) <- keep_cols

# labels (0 = germline, 1 = somatic)
y <- factor(df$sv_class, levels = c(1, 0), labels = c("somatic", "germline"))

tmp_env <- new.env()

# Load the file into the temporary environment
load("/data/beroukhim1/oumayma/no_bkpt_GaTSV/GaTS_BpA_031826.rda", envir = tmp_env)

# Extract the first (and likely only) object inside it and name it GaTSV
GaTSV <- tmp_env[[ls(tmp_env)[1]]]

# Optional: Clean up the temporary environment
rm(tmp_env)

# define results
results_list <- list()

# loop for iterating through seeds
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  cat("\nRunning replicate", i, "with seed", seeds[i], "...\n")
  
  # split dataset (2:1 train-test)
  trainIndex <- createDataPartition(y, p = 0.66, list = FALSE)
  X_train <- X[trainIndex, ]
  y_train <- y[trainIndex]
  X_test  <- X[-trainIndex, ]
  y_test  <- y[-trainIndex]
  
  # predicted probabilities on train set (Fixed model name to GaTSV)
  # predicted probabilities on train set
  prob_mat_train <- attr(predict(GaTSV, X_train, probability = TRUE), "probabilities")
  pos_col <- grep("somatic|1", colnames(prob_mat_train), value = TRUE, ignore.case = TRUE)
  if(length(pos_col) == 0) pos_col <- colnames(prob_mat_train)[2] else pos_col <- pos_col[1]
  probs_train <- as.numeric(prob_mat_train[, pos_col])  
  # get optimal cutoff
  cutoffs <- seq(0, 1, by = 0.001)
  best_cutoff <- 0
  best_score <- -Inf
  
  for (c in cutoffs) {
    preds <- factor(ifelse(probs_train > c, "somatic", "germline"), 
                    levels = c("somatic", "germline"))
    
    ppv <- ppv_vec(y_train, preds, event_level = "first")
    tpr <- sens_vec(y_train, preds, event_level = "first")
    
    if (is.na(ppv) || is.na(tpr)) next
    
    score <- ppv + tpr
    if (score > best_score) {
      best_score <- score
      best_cutoff <- c
    }
  }
  
  # evaluate on test set (Fixed model name to GaTSV)
  # evaluate on test set
  prob_mat_test <- attr(predict(GaTSV, X_test, probability = TRUE), "probabilities")
  probs_test <- as.numeric(prob_mat_test[, pos_col]) # Re-use the pos_col we found above!
  preds_test <- factor(ifelse(probs_test > best_cutoff, "somatic", "germline"), 
                       levels = c("somatic", "germline"))
  
  results_tbl <- tibble(truth = y_test, estimate = preds_test, prob = probs_test)
  
  # yardstick stats
  tpr_test  <- sens_vec(results_tbl$truth, results_tbl$estimate, event_level = "first")
  spec_test <- spec_vec(results_tbl$truth, results_tbl$estimate, event_level = "first")
  ppv_test  <- ppv_vec(results_tbl$truth, results_tbl$estimate, event_level = "first")
  auc_test  <- roc_auc_vec(results_tbl$truth, results_tbl$prob, event_level = "first")
  
  # store results
  results_list[[i]] <- tibble(
    replicate = i,
    seed = seeds[i],
    cutoff = best_cutoff,
    AUC = auc_test,
    PPV = ppv_test,
    TPR = tpr_test,
    Specificity = spec_test
  )
}

# save results
results_df <- bind_rows(results_list)
cat("\nFinal results:\n")
print(results_df)

cat("\nSummary statistics:\n")
print(summary(results_df))

# Save as an RDS file
saveRDS(results_df, "/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/svm_replicate_results.rds")